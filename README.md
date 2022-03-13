# Device-Free-Sensing-in-OFDM-Cellular-Network
The code is for the paper: [Q. Shi, L. Liu, S. Zhang, and S. Cui, "Device-free sensing in OFDM cellular network," to appear in IEEE J. Sel. Areas Commun..](https://arxiv.org/abs/2108.09177)
# Abstract
This paper considers device-free sensing in an orthogonal frequency division multiplexing (OFDM) cellular network to enable integrated sensing and communication (ISAC). A novel two-phase sensing framework is proposed to localize the passive targets that cannot transmit/receive reference signals to/from the base stations (BSs), where the ranges of the targets are estimated based on their reflected OFDM signals to the BSs in Phase I, and the location of each target is estimated based on its ranges to different BSs in Phase II. Specifically, in Phase I, we design a model-free range estimation approach by leveraging the OFDM channel estimation technique for determining the delay values of all the two-way BS-target-BS paths, which does not rely on any BS-target channel model. In Phase II, we reveal that ghost targets may be falsely detected in some cases as all the targets reflect the same signals to the BSs, which thus do not know how to match each estimated range with the right target. Interestingly, we show that the above data association issue is not a fundamental limitation for device-free sensing: under the ideal case of perfect range estimation in Phase I, the probability for ghost targets to exist is proved to be negligible when the targets are randomly located. Moreover, under the practical case of imperfect range estimation in Phase I, we propose an efficient algorithm for joint data association and target localization in Phase II. Numerical results show that our proposed two-phase framework can achieve very high accuracy in the localization of passive targets, which increases with the system bandwidth.
# Read Me
Model.m is the main function.\
coorest.m is used to emstimate target location by solving nonlinear sqaure problem through Gauss-Newton Method\
munkres.m is an open file [munkres](https://ww2.mathworks.cn/matlabcentral/fileexchange/20652-hungarian-algorithm-for-linear-assignment-problems-v2-3) in Matlab, which is used to solve two dimensional assignment problem through Hungarian Method.\
initial_est.m is the initial target estimation in the coorest.m.\
distsum.m is used to caculate the sum of distances from the estimated targets to all BSs.\
diffdist.m is used to calculate the mapping between the estimated targets and the true targets and the distances between them.

# Citation 
@article{shi2021device,\
  title={Device-Free Sensing in OFDM Cellular Network},\
  author={Shi, Qin and Liu, Liang and Zhang, Shuowen and Cui, Shuguang},\
  journal={arXiv preprint arXiv:2108.09177},\
  year={2021}\
}
# Note
The code is provided for the benefit of better understanding the results, and is not meant to be used in production.
