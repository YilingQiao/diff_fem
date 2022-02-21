import pydifem as pd
import torch
from api_diff import sim_all, sim_act
import time
import os
import argparse

parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('--inp', default='bridgesus', help='input json')
parser.add_argument('--out', default='bridgesus_demo', help='output dir')

args = parser.parse_args()

def get_loss(q, target):
    q = torch.reshape(q, [-1, 3])
    q = q.mean(0)
    print("simulate ans {}, target {}".format(q.data, target))
    loss = (q-target).norm()
    return loss 

def main():
    json_path = "/data/bridgesus/{}.json".format(args.inp)
    log_path = "out_test/{}".format(args.out)
    if not os.path.exists(log_path):
        os.mkdir(log_path)
    fo = open("{}/loss.txt".format(log_path), "w")
    
    num_sim_step = 10
    num_opt_step = 100

    env = pd.init_env(num_sim_step, json_path)
    init_x = pd.vector2double(env.mSoftWorld.mX)
    init_v = pd.vector2double(env.mSoftWorld.mV)
    init_jq = pd.vector2double(env.mSoftWorld.mJointQ)

    init_x = torch.tensor(init_x, dtype=torch.float32, requires_grad=True)
    init_v = torch.tensor(init_v, dtype=torch.float32, requires_grad=True)
    init_jq = torch.tensor(init_jq, dtype=torch.float32, requires_grad=True)
    init_tau = torch.tensor([1, 1]*668, dtype=torch.float32, requires_grad=True)
    t_scale = torch.tensor([2e5, 0.4]*668, dtype=torch.float32)
    optimizer = torch.optim.Adam([init_tau], lr=0.2)
    ini_mean = torch.reshape(init_x, [-1, 3]).mean(0)
    target = ini_mean + torch.tensor([0, -0.05, 0], dtype=torch.float32)

    print("init mean ", ini_mean)
    for steps in range(num_opt_step):
        optimizer.zero_grad()

        tau = init_tau
        tau = t_scale*torch.maximum(tau, tau*0+0.01)
        s = time.time()
        q, qd = sim_all(init_x, init_v, tau, env)
        
        print("forward time: {}".format(time.time() - s))
        loss = get_loss(q, target)
        print("---------------------------- iteration {:03d}: loss {}".format(steps, loss))
        s = time.time()
        loss.backward()
        print("backward time: {}".format(time.time() - s))
        optimizer.step()
        env.mPhase = 0
        fo.write( "{}\n".format(loss))
    fo.close()



if __name__ == "__main__":
  main()
