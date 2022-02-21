import torch
import pydifem as pd

class SimAct(torch.autograd.Function):
    @staticmethod
    def forward(ctx, q, qd, jq, jqd, tau, env, is_act):
        ctx.save_for_backward(q, qd, jq, jqd, tau)
        ctx.env = env
        ctx.phase = env.mPhase
        ctx.actionlevel = pd.vector2double(env.mDfobjs[0].GetActivationLevelsAggregate()) 
        ctx.is_act = is_act
        # print(ctx.actionlevel)
        ans = pd.forward_act(q.detach().numpy(), qd.detach().numpy(), 
            jq.detach().numpy(), jqd.detach().numpy(), tau.detach().numpy(), env, is_act)

        q, qd, jq, jqd = ans
        q = torch.tensor(q, dtype=torch.float32)
        qd = torch.tensor(qd, dtype=torch.float32)
        jq = torch.tensor(jq, dtype=torch.float32)
        jqd = torch.tensor(jqd, dtype=torch.float32)
        ans = q, qd, jq, jqd
        return q, qd, jq, jqd

    @staticmethod
    def backward(ctx, dldq, dldqd, dldjq, dldjqd):
        q, qd, jq, jqd, tau = ctx.saved_tensors
        ctx.env.mPhase = ctx.phase
        is_act = ctx.is_act 
        ctx.env.mDfobjs[0].SetActivationLevelsAggregate(
            pd.double2vector(ctx.actionlevel))

        # print(ctx.actionlevel)
        ans = pd.backward_act(q.detach().numpy(), qd.detach().numpy(), 
            jq.detach().numpy(), jqd.detach().numpy(), tau.detach().numpy(), dldq.detach().numpy(), 
            dldqd.detach().numpy(), dldjq.detach().numpy(), dldjqd.detach().numpy(), ctx.env, is_act)
        
        q, qd, jq, jqd, tau = ans
        q = torch.tensor(q, dtype=torch.float32)
        qd = torch.tensor(qd, dtype=torch.float32)
        jq = torch.tensor(jq, dtype=torch.float32)
        jqd = torch.tensor(jqd, dtype=torch.float32)
        tau = torch.tensor(tau, dtype=torch.float32)

        return q, qd, jq, jqd, tau, None, None

sim_act = SimAct.apply

class SimAll(torch.autograd.Function):
    @staticmethod
    def forward(ctx, q, qd, tau, env):
        ctx.env = env
        ans = pd.forward_all(q.detach().numpy(), qd.detach().numpy(), 
            tau.detach().numpy(), env)
        q, qd = ans
        ctx.save_for_backward(torch.tensor(q), torch.tensor(qd), tau)
        q = torch.tensor(q, dtype=torch.float32)
        qd = torch.tensor(qd, dtype=torch.float32)
        ans = q, qd
        return q, qd

    @staticmethod
    def backward(ctx, dldq, dldqd):
        q, qd, tau = ctx.saved_tensors
        ans = pd.backward_all(q.detach().numpy(), qd.detach().numpy(), 
            tau.detach().numpy(), dldq.detach().numpy(), dldqd.detach().numpy(), ctx.env)
        q, qd, tau = ans
        q = torch.tensor(q, dtype=torch.float32)
        qd = torch.tensor(qd, dtype=torch.float32)
        tau = torch.tensor(tau, dtype=torch.float32)
        return q, qd, tau, None

sim_all = SimAll.apply