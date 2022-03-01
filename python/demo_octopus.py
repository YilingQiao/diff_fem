import torch
import os
import time
from api_diff import sim_all, sim_act
import pydifem as pd
import argparse

parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('--inp', default='octopus', help='input json')
parser.add_argument('--out', default='octopus_demo', help='output dir')

args = parser.parse_args()

def get_loss(q, target):
    q = torch.reshape(q, [-1, 3])
    q = q.mean(0)
    print("simulate ans {}, target {}".format(q.data, target))
    loss = (q-target).norm()
    return loss 

def main():
    json_path = "/data/octopus/{}.json".format(args.inp)
    log_path = "out_test/{}".format(args.out)
    if not os.path.exists(log_path):
        os.mkdir(log_path)
    fo = open("{}/loss.txt".format(log_path), "w")

    num_sim_step = 400
    num_opt_step = 200
    simulation_Hz = 200

    env = pd.init_env(num_sim_step, json_path)
    init_x = pd.vector2double(env.mSoftWorld.mX)
    init_v = pd.vector2double(env.mSoftWorld.mV)
    init_jq = pd.vector2double(env.mSoftWorld.mJointQ)
    init_jqd = init_jq
    init_tau = [0] * 16 * 4
    print("init_jq: ", init_jq)

    init_x = torch.tensor(init_x, dtype=torch.float32, requires_grad=True)
    init_v = torch.tensor(init_v, dtype=torch.float32, requires_grad=True)
    init_jq = torch.tensor(init_jq, dtype=torch.float32, requires_grad=True)
    init_jqd = torch.tensor(init_jqd, dtype=torch.float32, requires_grad=True)
    init_tau = torch.tensor(init_tau*num_sim_step, dtype=torch.float32, requires_grad=True)
    taus = torch.reshape(init_tau, [num_sim_step, -1])
    print(taus.shape, " taus")
    optimizer = torch.optim.Adam([init_tau], lr=0.05)
    target = [0.4, 0.8, 0.4]
    target = torch.tensor(target, dtype=torch.float32, requires_grad=True)

    for steps in range(num_opt_step):
        optimizer.zero_grad()
        s = time.time()
        q, qd, jq, jqd = init_x, init_v, init_jq, init_jqd
        for i in range(num_sim_step):
            is_act = (i%6 == 0)
            q, qd, jq, jqd = sim_act(q, qd, jq, jqd, taus[i], env, is_act)
            
            env.SaveObj("{}/{:03d}_{:03d}".format(log_path, steps, i), 
                pd.double2vector(q.detach().numpy()))
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
