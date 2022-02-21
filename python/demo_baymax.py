import pydifem as pd
import torch
from api_diff import sim_all, sim_act
import time
import os
import numpy as np
import argparse
import sys

parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('--inp', default='baymax_skeleton', help='input json')
parser.add_argument('--out', default='baymax_demo', help='output dir')

args = parser.parse_args()

def get_loss(q, target, tau):
    q = torch.reshape(q, [-1, 3])
    q = q.mean(0)
    target_loss = (target[2]-q[2])**2
    smooth_loss =  args.w * (tau**2).mean()
    print("simulate ans {}, target loss {} smooth loss {}".format(
        q.data, target_loss, smooth_loss))
    loss = target_loss + smooth_loss
    return loss 

def main():
    json_path = "/data/baymax/{}.json".format(args.inp)
    log_path = "out_test/{}".format(args.out)
    if not os.path.exists(log_path):
        os.mkdir(log_path)
    
    num_sim_step = 100
    num_opt_step = 1
    simulationHz = 100

    env = pd.init_env(num_sim_step, json_path)
    init_x = pd.vector2double(env.mSoftWorld.mX)
    init_v = pd.vector2double(env.mSoftWorld.mV)
    init_jq = pd.vector2double(env.mSoftWorld.mJointQ)
    init_jqd = pd.vector2double(env.mSoftWorld.mJointQd)
    init_tau = pd.vector2double(env.mSoftWorld.mJointTorque)

    init_x = torch.tensor(init_x, dtype=torch.float32)
    init_v = torch.tensor(init_v, dtype=torch.float32)
    init_jq = torch.tensor(init_jq, dtype=torch.float32)
    init_jqd = torch.tensor(init_jqd, dtype=torch.float32)
    init_torque = torch.zeros([4], dtype=torch.float32)

    env.SaveObj("{}/000_000_init.obj".format(log_path), 
        pd.double2vector(init_x.detach().numpy()))

    for steps in range(num_opt_step):
        s = time.time()
        q, qd, jq, jqd = init_x, init_v, init_jq, init_jqd
        for i in range(num_sim_step):
            print("step: ", i)
            is_act = 1
            tau = init_torque
            tau[2] = 20
            tau[3] = -20
            if jq[10] < -0.66:
                tau[2] = 150
                tau[3] = -150
            if jq[14]>0.8 or jq[15]<-0.8:
                tau[2] = 0
                tau[3] = 0

            q, qd, jq, jqd = sim_act(q, qd, jq, jqd, tau, env, is_act)
          
            env.SaveObj("{}/{:03d}_{:03d}.obj".format(log_path, steps, i), 
                pd.double2vector(q.detach().numpy()))

if __name__ == "__main__":
  main()