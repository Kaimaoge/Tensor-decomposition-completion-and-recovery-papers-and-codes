# Collections of tensor reconstruction papers
Tensor reconstruction = Tensor completion and recovery

## Review of tensor decomposition
Kolda, Tamara G., and Brett W. Bader. "Tensor decompositions and applications." SIAM review 51.3 (2009): 455-500. [(Paper)](https://epubs.siam.org/doi/abs/10.1137/07070111X?journalCode=siread)

Cichocki, Andrzej, et al. Nonnegative matrix and tensor factorizations: applications to exploratory multi-way data analysis and blind source separation. John Wiley & Sons, 2009. [(Book)](http://www.academia.edu/download/46426061/Nonnegative_Matrix_and_Tensor_Factorizat20160612-12469-1usk837.pdf) 

### Big data and new decomposition algorithm
Cichocki, Andrzej. "Era of big data processing: A new approach via tensor networks and tensor decompositions." arXiv preprint arXiv:1403.2048 (2014).[(Paper)](https://arxiv.org/pdf/1403.2048.pdf) 

### Tensor decomposition and data fusion
Lahat, Dana, Tülay Adali, and Christian Jutten. "Multimodal data fusion: an overview of methods, challenges, and prospects." Proceedings of the IEEE 103.9 (2015): 1449-1477.[(Paper)](https://mdsoar.org/bitstream/handle/11603/10867/Lahat_Adali_Jutten_DataFusion_2015.pdf?sequence=1&isAllowed=y) 

### A good paper to understand tensor networks
Orús, Román. "A practical introduction to tensor networks: Matrix product states and projected entangled pair states." Annals of Physics 349 (2014): 117-158.[(Paper)](https://arxiv.org/pdf/1306.2164.pdf;)

## Tensor completion

### First definition of tensor nuclear norm, and transformation of tensor completion into a convex problem
Liu, Ji, et al. "Tensor completion for estimating missing values in visual data." IEEE transactions on pattern analysis and machine intelligence 35.1 (2012): 208-220.[(Paper)](https://repository.kaust.edu.sa/bitstream/handle/10754/562566/2012.PAMI.JiLiu.Tensor%20Completion.pdf?sequence=1) [(Code)](http://www.cs.rochester.edu/u/jliu/code/TensorCompletion.zip)

#### A parallel work
Gandy, Silvia, Benjamin Recht, and Isao Yamada. "Tensor completion and low-n-rank tensor recovery via convex optimization." Inverse Problems 27.2 (2011): 025010. [(Paper)](https://arxiv.org/pdf/1311.6182)

#### Trace norm of square matrix reshape of tensor (Theoritically more suitable for high dimensional tensor)
Mu, Cun, et al. "Square deal: Lower bounds and improved relaxations for tensor recovery." International conference on machine learning. 2014. [(Paper)](http://proceedings.mlr.press/v32/mu14.pdf)

### Tensor decomposition and completion

Acar, Evrim, et al. "Scalable tensor factorizations for incomplete data." Chemometrics and Intelligent Laboratory Systems 106.1 (2011): 41-56.[(Paper)](https://arxiv.org/pdf/1005.2197) [(Code)](https://www.tensortoolbox.org/cp_wopt_doc.html)- It has been a function of tensor toolbox

#### tensor singular value decomposition (t-SVD)
Zhang, Zemin, Gregory Ely, Shuchin Aeron, Ning Hao, and Misha Kilmer. "Novel methods for multilinear data completion and de-noising based on tensor-SVD." In Proceedings of the IEEE conference on computer vision and pattern recognition, pp. 3842-3849. 2014.[(Paper)](http://openaccess.thecvf.com/content_cvpr_2014/papers/Zhang_Novel_Methods_for_2014_CVPR_paper.pdf)

#### Simultaneous decomposition and completion (nuclear norm as a prior)
Chen, Yi-Lei, Chiou-Ting Hsu, and Hong-Yuan Mark Liao. "Simultaneous tensor decomposition and completion using factor priors." IEEE transactions on pattern analysis and machine intelligence 36.3 (2013): 577-591. [(Paper)](https://ir.nctu.edu.tw/bitstream/11536/23758/1/000331450100014.pdf) [(Code)](https://www.computer.org/csdl/journal/tp/2014/03/ttp2014030577/13rRUxASuNM)-supplementary material

#### Bayesian CP factorization for incomplete tensor
Zhao, Qibin, Liqing Zhang, and Andrzej Cichocki. "Bayesian CP factorization of incomplete tensors with automatic rank determination." IEEE transactions on pattern analysis and machine intelligence 37.9 (2015): 1751-1763. [(Paper)](https://arxiv.org/pdf/1401.6497) [(Code)](https://github.com/qbzhao/BCPF)

#### Smooth prior for CP factorization (Achieves the best performance in my comparison experiments)
Yokota, Tatsuya, Qibin Zhao, and Andrzej Cichocki. "Smooth PARAFAC decomposition for tensor completion." IEEE Transactions on Signal Processing 64.20 (2016): 5423-5436.[(Paper)](https://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=7502115) [(Code)](https://drive.google.com/file/d/15xk67wYZ9GI2Kn93g_aaA4CsmgnqGlHE/view)

#### A tensor completion framework fused multiple priors
Wu, Y., Tan, H., Li, Y., Zhang, J., & Chen, X. (2018). A fused CP factorization method for incomplete tensors. IEEE transactions on neural networks and learning systems, 30(3), 751-764. [(Nonnegative factorization Code)](https://github.com/Kaimaoge/Collections-of-tensor-reconstruction-papers/tree/master/NonnegativeFCP)

#### Beyond low-rank, mixture gaussian for tensor modeling
Zhang, L., Wei, W., Shi, Q., Shen, C., Hengel, A. V. D., & Zhang, Y. (2017). Beyond low rank: A data-adaptive tensor completion method. arXiv preprint arXiv:1708.01008.[(Paper)](https://arxiv.org/pdf/1708.01008.pdf)


## New decomposition

### Tensor ring decomposition (Defined upon trace operation)
Zhao, Qibin, et al. "Tensor ring decomposition." arXiv preprint arXiv:1606.05535 (2016).[(Paper)](https://arxiv.org/pdf/1606.05535.pdf)[(Matlab Code)](https://github.com/oscarmickelin/tensor-ring-decomposition)

### Tensor train decomposition (analogy to matrix product state in quantum physics)
Oseledets, Ivan V. "Tensor-train decomposition." SIAM Journal on Scientific Computing 33.5 (2011): 2295-2317. [(Paper)](https://www.researchgate.net/profile/Ivan_Oseledets2/publication/220412263_Tensor-Train_Decomposition/links/5bbfb5c5299bf1004c5a56e3/Tensor-Train-Decomposition.pdf) [(Matlab code)](https://github.com/oseledets/TT-Toolbox) [(Tensorflow code)](https://github.com/Bihaqo/t3f)


### To be continued

