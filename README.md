
# Differentiable Simulation of Soft Multi-body Systems

*Yi-Ling Qiao, Junbang Liang, Vladlen Koltun, Ming C. Lin*

[[Paper]](http://vladlen.info/publications/differentiable-simulation-soft-multi-body-systems/) [[Code]](https://github.com/YilingQiao/diff_fem) [[Video]](https://youtu.be/TPgFM5WxzaU)
## Updates
1. Link to a demo video.
2. cmakelist and setup files.
3. demos for inverse problems and forward simulation.

## TODO
1. More demos for the control problems.
2. More documents for the config files.
3. Utils for make tet meshes.
4. More readme documentation.

## Setup
1. Create a conda virtual environment and activate it.
```bash
conda create -n difem python=3.8 -y
conda activate difem
```

2. Download and build the project.
```bash
git clone git@github.com:YilingQiao/diff_fem.git
cd diff_fem
git submodule init
git submodule update
sudo apt-get install ninja-build cppad libcgal-dev
python setup.py install
```
3. Run the examples
## Examples
### Inverse problem
1. Suspension bridge (Fig. 3a in the paper)
```bash
python python/demo_sus.py
```
2. Arch bridge (Fig. 3c in the paper)
```bash
python python/demo_br.py
```

For the above 2 experiments, the output meshes are stored in `out/`
### Control Problems
1. Drone (Fig. 4a in the paper)
TODO

2. Octopus (Fig. 4b in the paper)
```bash
python python/demo_octopus.py
```

3. Fish (Fig. 4c in the paper)
TODO

For the above experiments, the output meshes are stored in `out_test/`
### Forward Simulation
Note that our simulator can be used for pure forward simulation. In this case, we replace the autodiff scalar (cppad) by C++ double and can run much faster (more than 5x).

To make this change, we first uncomment 
```cpp
// #define FORWARD_ONLY
```
in `python/pydifem.cc` and then run 
```bash
python setup.py install
```

1. Baymax (Fig. 1 in the paper)
```bash
python python/demo_baymax.py
```
2. Clothball (Fig. 2 in the paper)
```bash
python python/demo_cloth_ball.py
```

For the above experiments, the output meshes are stored in `out_test/`
## Our Related Repos
Differentiable Soft Body Dynamics (this repo) [Code](https://github.com/YilingQiao/diff_fem) [Paper](http://vladlen.info/publications/differentiable-simulation-soft-multi-body-systems/)
*Differentiable Simulation of Soft Multi-body Systems. Yi-Ling Qiao, Junbang Liang, Vladlen Koltun, Ming C. Lin. (Neurips 2021)*

Differentiable Articulated Body Dynamics [Code](https://github.com/YilingQiao/diffarticulated) [Paper](https://arxiv.org/abs/2109.07719)
*Efficient Differentiable Simulation of Articulated Bodies. Yi-Ling Qiao, Junbang Liang, Vladlen Koltun, Ming C. Lin. (ICML 2021)*

Differentiable Dynamics for Rigid Body and Cloth Coupling [Code](https://github.com/YilingQiao/diffsim) [Paper](https://arxiv.org/abs/2007.02168)
*Scalable Differentiable Physics for Learning and Control. Yi-Ling Qiao, Junbang Liang, Vladlen Koltun, Ming C. Lin. (ICML 2020)*

Differentiable Cloth Dynamics [Code](https://github.com/williamljb/DifferentiableCloth) [Paper](https://www.cs.umd.edu/~liangjb/docs/NIPS2019.pdf)
*Differentiable Cloth Simulation for Inverse Problems. Junbang Liang, Ming C. Lin, Vladlen Koltun. (NeurIPS 2019)*

## Bibtex
```
@inproceedings{Qiao2021Differentiable,
author  = {Qiao, Yi-Ling and Liang, Junbang and Koltun, Vladlen and Lin, Ming C.},
title  = {Differentiable Simulation of Soft Multi-body Systems},
booktitle = {Conference on Neural Information Processing Systems (NeurIPS)},
year  = {2021},
}
```
