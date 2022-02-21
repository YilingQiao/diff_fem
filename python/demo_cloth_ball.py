import pydifem as pd
import torch
from api_diff import sim_all, sim_act
import time
import os
import numpy as np
import argparse
import sys

parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('--inp', default='cloth_ball', help='input json')
parser.add_argument('--out', default='cloth_ball_demo', help='output dir')

args = parser.parse_args()

def main():
    json_path = "/data/cloth_ball/{}.json".format(args.inp)
    log_path = "out_test/{}".format(args.out)
    if not os.path.exists(log_path):
        os.mkdir(log_path)
    
    num_sim_step = 100
    num_opt_step = 1

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

    env.SaveObj("{}/000_000_init".format(log_path), 
        pd.double2vector(init_x.detach().numpy()))

    for steps in range(num_opt_step):
        q, qd, jq, jqd = init_x, init_v, init_jq, init_jqd
        for i in range(num_sim_step):
            print("step: ", i)
            is_act = 0
            tau = jq
            q, qd, jq, jqd = sim_act(q, qd, jq, jqd, tau, env, is_act)
          
            env.SaveObj("{}/{:03d}_{:03d}".format(log_path, steps, i), 
                pd.double2vector(q.detach().numpy()))

if __name__ == "__main__":
  main()