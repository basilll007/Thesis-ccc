## Phase 2 baselines (grouped 5-fold, pooled held-out R2)

| receiver   | target            | model          | R2 [95% tumour CI]     | gain vs own covariates [CI]   | config                                          |
|:-----------|:------------------|:---------------|:-----------------------|:------------------------------|:------------------------------------------------|
| myeloid    | C1QA              | mean/none      | 0.000 [0.000, 0.000]   | -0.180 [-0.227, -0.101]       | train_mean                                      |
| myeloid    | C1QB              | mean/none      | 0.000 [0.000, 0.000]   | -0.255 [-0.304, -0.175]       | train_mean                                      |
| myeloid    | C1QC              | mean/none      | 0.000 [0.000, 0.000]   | -0.237 [-0.302, -0.139]       | train_mean                                      |
| myeloid    | CD163             | mean/none      | 0.000 [0.000, 0.000]   | -0.186 [-0.237, -0.106]       | train_mean                                      |
| myeloid    | CD9               | mean/none      | 0.000 [0.000, 0.000]   | -0.106 [-0.130, -0.070]       | train_mean                                      |
| myeloid    | FABP5             | mean/none      | 0.000 [0.000, 0.000]   | -0.118 [-0.142, -0.083]       | train_mean                                      |
| myeloid    | GPNMB             | mean/none      | 0.000 [0.000, 0.000]   | -0.193 [-0.224, -0.146]       | train_mean                                      |
| myeloid    | LGALS3            | mean/none      | 0.000 [0.000, 0.000]   | -0.136 [-0.160, -0.109]       | train_mean                                      |
| myeloid    | MRC1              | mean/none      | 0.000 [0.000, 0.000]   | -0.094 [-0.118, -0.062]       | train_mean                                      |
| myeloid    | SELENOP           | mean/none      | 0.000 [0.000, 0.000]   | -0.229 [-0.273, -0.162]       | train_mean                                      |
| myeloid    | SPP1              | mean/none      | 0.000 [0.000, 0.000]   | -0.613 [-0.707, -0.478]       | train_mean                                      |
| myeloid    | genes_pooled      | mean/none      | 0.000 [0.000, 0.000]   | -0.265 [-0.320, -0.191]       | train_mean                                      |
| myeloid    | myeloid:C1Q_TAM   | mean/none      | 0.000 [0.000, 0.000]   | -0.299 [-0.382, -0.188]       | train_mean                                      |
| myeloid    | myeloid:SPP1_TAM  | mean/none      | 0.000 [0.000, 0.000]   | -0.335 [-0.398, -0.254]       | train_mean                                      |
| myeloid    | C1QA              | mlp/genes      | 0.073 [-0.026, 0.145]  | -0.107 [-0.191, -0.041]       | lr=0.001,seed=0;lr=0.001,seed=1;lr=0.001,seed=2 |
| myeloid    | C1QB              | mlp/genes      | 0.125 [-0.005, 0.230]  | -0.131 [-0.238, -0.037]       | lr=0.001,seed=0;lr=0.001,seed=1;lr=0.001,seed=2 |
| myeloid    | C1QC              | mlp/genes      | 0.108 [-0.020, 0.210]  | -0.129 [-0.232, -0.051]       | lr=0.001,seed=0;lr=0.001,seed=1;lr=0.001,seed=2 |
| myeloid    | CD163             | mlp/genes      | 0.052 [-0.053, 0.130]  | -0.134 [-0.209, -0.073]       | lr=0.001,seed=0;lr=0.001,seed=1;lr=0.001,seed=2 |
| myeloid    | CD9               | mlp/genes      | 0.079 [0.056, 0.095]   | -0.027 [-0.039, -0.006]       | lr=0.001,seed=0;lr=0.001,seed=1;lr=0.001,seed=2 |
| myeloid    | FABP5             | mlp/genes      | 0.100 [0.059, 0.133]   | -0.018 [-0.033, -0.003]       | lr=0.001,seed=0;lr=0.001,seed=1;lr=0.001,seed=2 |
| myeloid    | GPNMB             | mlp/genes      | 0.134 [0.067, 0.172]   | -0.058 [-0.104, -0.032]       | lr=0.001,seed=0;lr=0.001,seed=1;lr=0.001,seed=2 |
| myeloid    | LGALS3            | mlp/genes      | 0.084 [0.054, 0.115]   | -0.052 [-0.079, -0.029]       | lr=0.001,seed=0;lr=0.001,seed=1;lr=0.001,seed=2 |
| myeloid    | MRC1              | mlp/genes      | 0.016 [-0.013, 0.056]  | -0.078 [-0.121, -0.031]       | lr=0.001,seed=0;lr=0.001,seed=1;lr=0.001,seed=2 |
| myeloid    | SELENOP           | mlp/genes      | 0.204 [0.136, 0.261]   | -0.025 [-0.061, 0.017]        | lr=0.001,seed=0;lr=0.001,seed=1;lr=0.001,seed=2 |
| myeloid    | SPP1              | mlp/genes      | 0.543 [0.333, 0.672]   | -0.070 [-0.156, -0.027]       | lr=0.001,seed=0;lr=0.001,seed=1;lr=0.001,seed=2 |
| myeloid    | genes_pooled      | mlp/genes      | 0.183 [0.087, 0.265]   | -0.082 [-0.135, -0.045]       | lr=0.001,seed=0;lr=0.001,seed=1;lr=0.001,seed=2 |
| myeloid    | myeloid:C1Q_TAM   | mlp/genes      | 0.187 [0.008, 0.316]   | -0.113 [-0.225, -0.045]       | lr=0.001,seed=0;lr=0.001,seed=1;lr=0.001,seed=2 |
| myeloid    | myeloid:SPP1_TAM  | mlp/genes      | 0.334 [0.234, 0.405]   | -0.000 [-0.028, 0.017]        | lr=0.001,seed=0;lr=0.001,seed=1;lr=0.001,seed=2 |
| myeloid    | C1QA              | ridge/genes    | 0.161 [0.091, 0.204]   | -0.019 [-0.033, -0.005]       | alpha=1000.0                                    |
| myeloid    | C1QB              | ridge/genes    | 0.251 [0.176, 0.301]   | -0.004 [-0.019, 0.012]        | alpha=1000.0                                    |
| myeloid    | C1QC              | ridge/genes    | 0.235 [0.136, 0.309]   | -0.003 [-0.020, 0.012]        | alpha=1000.0                                    |
| myeloid    | CD163             | ridge/genes    | 0.157 [0.091, 0.204]   | -0.029 [-0.048, -0.006]       | alpha=1000.0                                    |
| myeloid    | CD9               | ridge/genes    | 0.087 [0.041, 0.119]   | -0.019 [-0.034, -0.008]       | alpha=1000.0                                    |
| myeloid    | FABP5             | ridge/genes    | 0.094 [0.043, 0.130]   | -0.024 [-0.042, -0.011]       | alpha=1000.0                                    |
| myeloid    | GPNMB             | ridge/genes    | 0.191 [0.135, 0.232]   | -0.002 [-0.028, 0.016]        | alpha=1000.0                                    |
| myeloid    | LGALS3            | ridge/genes    | 0.106 [0.077, 0.133]   | -0.030 [-0.036, -0.023]       | alpha=1000.0                                    |
| myeloid    | MRC1              | ridge/genes    | 0.056 [0.038, 0.073]   | -0.038 [-0.056, -0.018]       | alpha=1000.0                                    |
| myeloid    | SELENOP           | ridge/genes    | 0.235 [0.174, 0.279]   | 0.006 [-0.010, 0.031]         | alpha=1000.0                                    |
| myeloid    | SPP1              | ridge/genes    | 0.617 [0.497, 0.700]   | 0.004 [-0.011, 0.038]         | alpha=1000.0                                    |
| myeloid    | genes_pooled      | ridge/genes    | 0.253 [0.185, 0.309]   | -0.012 [-0.023, 0.002]        | alpha=1000.0                                    |
| myeloid    | myeloid:C1Q_TAM   | ridge/genes    | 0.311 [0.204, 0.406]   | 0.012 [-0.014, 0.035]         | alpha=1000.0                                    |
| myeloid    | myeloid:SPP1_TAM  | ridge/genes    | 0.349 [0.265, 0.406]   | 0.014 [-0.003, 0.042]         | alpha=1000.0                                    |
| myeloid    | C1QA              | ridge/own      | 0.180 [0.101, 0.227]   | 0.000 [0.000, 0.000]          | alpha=1000.0                                    |
| myeloid    | C1QB              | ridge/own      | 0.255 [0.175, 0.304]   | 0.000 [0.000, 0.000]          | alpha=1000.0                                    |
| myeloid    | C1QC              | ridge/own      | 0.237 [0.139, 0.302]   | 0.000 [0.000, 0.000]          | alpha=1000.0                                    |
| myeloid    | CD163             | ridge/own      | 0.186 [0.106, 0.237]   | 0.000 [0.000, 0.000]          | alpha=1000.0                                    |
| myeloid    | CD9               | ridge/own      | 0.106 [0.070, 0.130]   | 0.000 [0.000, 0.000]          | alpha=1000.0                                    |
| myeloid    | FABP5             | ridge/own      | 0.118 [0.083, 0.142]   | 0.000 [0.000, 0.000]          | alpha=1000.0                                    |
| myeloid    | GPNMB             | ridge/own      | 0.193 [0.146, 0.224]   | 0.000 [0.000, 0.000]          | alpha=1000.0                                    |
| myeloid    | LGALS3            | ridge/own      | 0.136 [0.109, 0.160]   | 0.000 [0.000, 0.000]          | alpha=1000.0                                    |
| myeloid    | MRC1              | ridge/own      | 0.094 [0.062, 0.118]   | 0.000 [0.000, 0.000]          | alpha=1000.0                                    |
| myeloid    | SELENOP           | ridge/own      | 0.229 [0.162, 0.273]   | 0.000 [0.000, 0.000]          | alpha=1000.0                                    |
| myeloid    | SPP1              | ridge/own      | 0.613 [0.478, 0.707]   | 0.000 [0.000, 0.000]          | alpha=1000.0                                    |
| myeloid    | genes_pooled      | ridge/own      | 0.265 [0.191, 0.320]   | 0.000 [0.000, 0.000]          | alpha=1000.0                                    |
| myeloid    | myeloid:C1Q_TAM   | ridge/own      | 0.299 [0.188, 0.382]   | 0.000 [0.000, 0.000]          | alpha=1000.0                                    |
| myeloid    | myeloid:SPP1_TAM  | ridge/own      | 0.335 [0.254, 0.398]   | 0.000 [0.000, 0.000]          | alpha=1000.0                                    |
| myeloid    | C1QA              | ridge/programs | 0.197 [0.116, 0.244]   | 0.017 [0.011, 0.021]          | alpha=0.1                                       |
| myeloid    | C1QB              | ridge/programs | 0.261 [0.188, 0.306]   | 0.005 [-0.000, 0.015]         | alpha=0.1                                       |
| myeloid    | C1QC              | ridge/programs | 0.251 [0.151, 0.320]   | 0.013 [0.005, 0.020]          | alpha=0.1                                       |
| myeloid    | CD163             | ridge/programs | 0.190 [0.104, 0.247]   | 0.004 [-0.005, 0.011]         | alpha=0.1                                       |
| myeloid    | CD9               | ridge/programs | 0.112 [0.073, 0.137]   | 0.006 [0.002, 0.009]          | alpha=0.1                                       |
| myeloid    | FABP5             | ridge/programs | 0.124 [0.083, 0.153]   | 0.006 [-0.001, 0.013]         | alpha=0.1                                       |
| myeloid    | GPNMB             | ridge/programs | 0.212 [0.159, 0.251]   | 0.019 [0.009, 0.029]          | alpha=0.1                                       |
| myeloid    | LGALS3            | ridge/programs | 0.142 [0.114, 0.167]   | 0.006 [0.001, 0.011]          | alpha=0.1                                       |
| myeloid    | MRC1              | ridge/programs | 0.095 [0.064, 0.119]   | 0.001 [-0.003, 0.006]         | alpha=0.1                                       |
| myeloid    | SELENOP           | ridge/programs | 0.245 [0.178, 0.288]   | 0.016 [0.009, 0.025]          | alpha=0.1                                       |
| myeloid    | SPP1              | ridge/programs | 0.633 [0.498, 0.728]   | 0.020 [0.011, 0.031]          | alpha=0.1                                       |
| myeloid    | genes_pooled      | ridge/programs | 0.277 [0.201, 0.332]   | 0.012 [0.008, 0.016]          | alpha=0.1                                       |
| myeloid    | myeloid:C1Q_TAM   | ridge/programs | 0.333 [0.223, 0.417]   | 0.033 [0.023, 0.043]          | alpha=0.1                                       |
| myeloid    | myeloid:SPP1_TAM  | ridge/programs | 0.383 [0.299, 0.443]   | 0.048 [0.037, 0.057]          | alpha=0.1                                       |
| myeloid    | C1QA              | xgb/genes      | 0.000 [0.000, 0.000]   | -0.180 [-0.227, -0.101]       | lr=0.05,leaves=15                               |
| myeloid    | C1QB              | xgb/genes      | 0.000 [0.000, 0.000]   | -0.255 [-0.304, -0.175]       | lr=0.05,leaves=15                               |
| myeloid    | C1QC              | xgb/genes      | 0.000 [0.000, 0.000]   | -0.237 [-0.302, -0.139]       | lr=0.05,leaves=15                               |
| myeloid    | CD163             | xgb/genes      | 0.000 [0.000, 0.000]   | -0.186 [-0.237, -0.106]       | lr=0.05,leaves=15                               |
| myeloid    | CD9               | xgb/genes      | 0.000 [0.000, 0.000]   | -0.106 [-0.130, -0.070]       | lr=0.05,leaves=15                               |
| myeloid    | FABP5             | xgb/genes      | 0.000 [0.000, 0.000]   | -0.118 [-0.142, -0.083]       | lr=0.05,leaves=15                               |
| myeloid    | GPNMB             | xgb/genes      | 0.000 [0.000, 0.000]   | -0.193 [-0.224, -0.146]       | lr=0.05,leaves=15                               |
| myeloid    | LGALS3            | xgb/genes      | 0.000 [0.000, 0.000]   | -0.136 [-0.160, -0.109]       | lr=0.05,leaves=15                               |
| myeloid    | MRC1              | xgb/genes      | 0.000 [0.000, 0.000]   | -0.094 [-0.118, -0.062]       | lr=0.05,leaves=15                               |
| myeloid    | SELENOP           | xgb/genes      | 0.000 [0.000, 0.000]   | -0.229 [-0.273, -0.162]       | lr=0.05,leaves=15                               |
| myeloid    | SPP1              | xgb/genes      | 0.000 [0.000, 0.000]   | -0.613 [-0.707, -0.478]       | lr=0.05,leaves=15                               |
| myeloid    | genes_pooled      | xgb/genes      | 0.000 [0.000, 0.000]   | -0.265 [-0.320, -0.191]       | lr=0.05,leaves=15                               |
| myeloid    | myeloid:C1Q_TAM   | xgb/genes      | 0.334 [0.222, 0.424]   | 0.034 [0.018, 0.048]          | lr=0.05,leaves=15                               |
| myeloid    | myeloid:SPP1_TAM  | xgb/genes      | 0.385 [0.312, 0.436]   | 0.050 [0.035, 0.070]          | lr=0.05,leaves=15                               |
| myeloid    | C1QA              | xgb/programs   | 0.000 [0.000, 0.000]   | -0.180 [-0.227, -0.101]       | lr=0.1,leaves=15                                |
| myeloid    | C1QB              | xgb/programs   | 0.000 [0.000, 0.000]   | -0.255 [-0.304, -0.175]       | lr=0.1,leaves=15                                |
| myeloid    | C1QC              | xgb/programs   | 0.000 [0.000, 0.000]   | -0.237 [-0.302, -0.139]       | lr=0.1,leaves=15                                |
| myeloid    | CD163             | xgb/programs   | 0.000 [0.000, 0.000]   | -0.186 [-0.237, -0.106]       | lr=0.1,leaves=15                                |
| myeloid    | CD9               | xgb/programs   | 0.000 [0.000, 0.000]   | -0.106 [-0.130, -0.070]       | lr=0.1,leaves=15                                |
| myeloid    | FABP5             | xgb/programs   | 0.000 [0.000, 0.000]   | -0.118 [-0.142, -0.083]       | lr=0.1,leaves=15                                |
| myeloid    | GPNMB             | xgb/programs   | 0.000 [0.000, 0.000]   | -0.193 [-0.224, -0.146]       | lr=0.1,leaves=15                                |
| myeloid    | LGALS3            | xgb/programs   | 0.000 [0.000, 0.000]   | -0.136 [-0.160, -0.109]       | lr=0.1,leaves=15                                |
| myeloid    | MRC1              | xgb/programs   | 0.000 [0.000, 0.000]   | -0.094 [-0.118, -0.062]       | lr=0.1,leaves=15                                |
| myeloid    | SELENOP           | xgb/programs   | 0.000 [0.000, 0.000]   | -0.229 [-0.273, -0.162]       | lr=0.1,leaves=15                                |
| myeloid    | SPP1              | xgb/programs   | 0.000 [0.000, 0.000]   | -0.613 [-0.707, -0.478]       | lr=0.1,leaves=15                                |
| myeloid    | genes_pooled      | xgb/programs   | 0.000 [0.000, 0.000]   | -0.265 [-0.320, -0.191]       | lr=0.1,leaves=15                                |
| myeloid    | myeloid:C1Q_TAM   | xgb/programs   | 0.331 [0.222, 0.420]   | 0.032 [0.018, 0.044]          | lr=0.1,leaves=15                                |
| myeloid    | myeloid:SPP1_TAM  | xgb/programs   | 0.387 [0.307, 0.446]   | 0.053 [0.043, 0.064]          | lr=0.1,leaves=15                                |
| tumour     | AGR2              | mean/none      | 0.000 [0.000, 0.000]   | -0.451 [-0.520, -0.362]       | train_mean                                      |
| tumour     | AREG              | mean/none      | 0.000 [0.000, 0.000]   | -0.157 [-0.263, -0.059]       | train_mean                                      |
| tumour     | CEACAM1           | mean/none      | 0.000 [0.000, 0.000]   | -0.132 [-0.173, -0.095]       | train_mean                                      |
| tumour     | CEACAM6           | mean/none      | 0.000 [0.000, 0.000]   | -0.627 [-0.722, -0.480]       | train_mean                                      |
| tumour     | ITGA3             | mean/none      | 0.000 [0.000, 0.000]   | -0.238 [-0.295, -0.201]       | train_mean                                      |
| tumour     | KRT14             | mean/none      | 0.000 [0.000, 0.000]   | -0.307 [-0.519, -0.117]       | train_mean                                      |
| tumour     | KRT15             | mean/none      | 0.000 [0.000, 0.000]   | -0.089 [-0.109, -0.072]       | train_mean                                      |
| tumour     | KRT17             | mean/none      | 0.000 [0.000, 0.000]   | -0.428 [-0.548, -0.328]       | train_mean                                      |
| tumour     | KRT6A             | mean/none      | 0.000 [0.000, 0.000]   | -0.574 [-0.781, -0.290]       | train_mean                                      |
| tumour     | LCN2              | mean/none      | 0.000 [0.000, 0.000]   | -0.540 [-0.622, -0.447]       | train_mean                                      |
| tumour     | MET               | mean/none      | 0.000 [0.000, 0.000]   | -0.086 [-0.106, -0.065]       | train_mean                                      |
| tumour     | OLFM4             | mean/none      | 0.000 [0.000, 0.000]   | -0.604 [-0.736, -0.387]       | train_mean                                      |
| tumour     | PIGR              | mean/none      | 0.000 [0.000, 0.000]   | -0.535 [-0.628, -0.411]       | train_mean                                      |
| tumour     | S100A2            | mean/none      | 0.000 [0.000, 0.000]   | -0.453 [-0.648, -0.197]       | train_mean                                      |
| tumour     | S100P             | mean/none      | 0.000 [0.000, 0.000]   | -0.466 [-0.568, -0.358]       | train_mean                                      |
| tumour     | SERPINB5          | mean/none      | 0.000 [0.000, 0.000]   | -0.147 [-0.255, -0.083]       | train_mean                                      |
| tumour     | SLC2A1            | mean/none      | 0.000 [0.000, 0.000]   | -0.217 [-0.269, -0.167]       | train_mean                                      |
| tumour     | SPINK1            | mean/none      | 0.000 [0.000, 0.000]   | -0.579 [-0.656, -0.466]       | train_mean                                      |
| tumour     | genes_pooled      | mean/none      | 0.000 [0.000, 0.000]   | -0.425 [-0.514, -0.352]       | train_mean                                      |
| tumour     | tumour:basal_like | mean/none      | 0.000 [0.000, 0.000]   | -0.397 [-0.512, -0.289]       | train_mean                                      |
| tumour     | tumour:classical  | mean/none      | 0.000 [0.000, 0.000]   | -0.618 [-0.711, -0.502]       | train_mean                                      |
| tumour     | AGR2              | mlp/genes      | 0.290 [0.200, 0.409]   | -0.160 [-0.240, -0.060]       | lr=0.001,seed=0;lr=0.001,seed=1;lr=0.001,seed=2 |
| tumour     | AREG              | mlp/genes      | 0.082 [0.012, 0.136]   | -0.076 [-0.128, -0.032]       | lr=0.001,seed=0;lr=0.001,seed=1;lr=0.001,seed=2 |
| tumour     | CEACAM1           | mlp/genes      | 0.068 [0.027, 0.107]   | -0.065 [-0.089, -0.043]       | lr=0.001,seed=0;lr=0.001,seed=1;lr=0.001,seed=2 |
| tumour     | CEACAM6           | mlp/genes      | 0.405 [0.235, 0.562]   | -0.221 [-0.382, -0.092]       | lr=0.001,seed=0;lr=0.001,seed=1;lr=0.001,seed=2 |
| tumour     | ITGA3             | mlp/genes      | 0.048 [-0.027, 0.129]  | -0.191 [-0.241, -0.145]       | lr=0.001,seed=0;lr=0.001,seed=1;lr=0.001,seed=2 |
| tumour     | KRT14             | mlp/genes      | 0.087 [-0.030, 0.180]  | -0.219 [-0.358, -0.091]       | lr=0.001,seed=0;lr=0.001,seed=1;lr=0.001,seed=2 |
| tumour     | KRT15             | mlp/genes      | -0.026 [-0.063, 0.018] | -0.116 [-0.147, -0.081]       | lr=0.001,seed=0;lr=0.001,seed=1;lr=0.001,seed=2 |
| tumour     | KRT17             | mlp/genes      | 0.016 [-0.142, 0.158]  | -0.412 [-0.573, -0.261]       | lr=0.001,seed=0;lr=0.001,seed=1;lr=0.001,seed=2 |
| tumour     | KRT6A             | mlp/genes      | 0.488 [0.200, 0.696]   | -0.085 [-0.111, -0.057]       | lr=0.001,seed=0;lr=0.001,seed=1;lr=0.001,seed=2 |
| tumour     | LCN2              | mlp/genes      | 0.275 [0.098, 0.453]   | -0.265 [-0.371, -0.149]       | lr=0.001,seed=0;lr=0.001,seed=1;lr=0.001,seed=2 |
| tumour     | MET               | mlp/genes      | 0.019 [-0.011, 0.047]  | -0.067 [-0.099, -0.039]       | lr=0.001,seed=0;lr=0.001,seed=1;lr=0.001,seed=2 |
| tumour     | OLFM4             | mlp/genes      | 0.512 [0.297, 0.657]   | -0.092 [-0.136, -0.050]       | lr=0.001,seed=0;lr=0.001,seed=1;lr=0.001,seed=2 |
| tumour     | PIGR              | mlp/genes      | 0.467 [0.335, 0.572]   | -0.068 [-0.092, -0.045]       | lr=0.001,seed=0;lr=0.001,seed=1;lr=0.001,seed=2 |
| tumour     | S100A2            | mlp/genes      | 0.313 [0.097, 0.517]   | -0.139 [-0.219, -0.042]       | lr=0.001,seed=0;lr=0.001,seed=1;lr=0.001,seed=2 |
| tumour     | S100P             | mlp/genes      | 0.210 [0.076, 0.367]   | -0.256 [-0.341, -0.157]       | lr=0.001,seed=0;lr=0.001,seed=1;lr=0.001,seed=2 |
| tumour     | SERPINB5          | mlp/genes      | 0.065 [-0.015, 0.176]  | -0.082 [-0.115, -0.046]       | lr=0.001,seed=0;lr=0.001,seed=1;lr=0.001,seed=2 |
| tumour     | SLC2A1            | mlp/genes      | 0.087 [0.024, 0.145]   | -0.131 [-0.183, -0.082]       | lr=0.001,seed=0;lr=0.001,seed=1;lr=0.001,seed=2 |
| tumour     | SPINK1            | mlp/genes      | 0.374 [0.273, 0.501]   | -0.205 [-0.292, -0.093]       | lr=0.001,seed=0;lr=0.001,seed=1;lr=0.001,seed=2 |
| tumour     | genes_pooled      | mlp/genes      | 0.252 [0.177, 0.353]   | -0.173 [-0.209, -0.133]       | lr=0.001,seed=0;lr=0.001,seed=1;lr=0.001,seed=2 |
| tumour     | tumour:basal_like | mlp/genes      | 0.265 [0.126, 0.406]   | -0.133 [-0.222, -0.071]       | lr=0.001,seed=0;lr=0.001,seed=1;lr=0.001,seed=2 |
| tumour     | tumour:classical  | mlp/genes      | 0.547 [0.416, 0.644]   | -0.071 [-0.132, -0.021]       | lr=0.001,seed=0;lr=0.001,seed=1;lr=0.001,seed=2 |
| tumour     | AGR2              | ridge/genes    | 0.441 [0.347, 0.512]   | -0.010 [-0.033, 0.009]        | alpha=1000.0                                    |
| tumour     | AREG              | ridge/genes    | 0.085 [-0.034, 0.188]  | -0.073 [-0.109, -0.045]       | alpha=1000.0                                    |
| tumour     | CEACAM1           | ridge/genes    | 0.091 [0.029, 0.152]   | -0.041 [-0.070, -0.019]       | alpha=1000.0                                    |
| tumour     | CEACAM6           | ridge/genes    | 0.646 [0.524, 0.728]   | 0.020 [-0.011, 0.054]         | alpha=1000.0                                    |
| tumour     | ITGA3             | ridge/genes    | 0.201 [0.176, 0.227]   | -0.038 [-0.077, -0.012]       | alpha=1000.0                                    |
| tumour     | KRT14             | ridge/genes    | 0.255 [0.088, 0.428]   | -0.052 [-0.093, -0.014]       | alpha=1000.0                                    |
| tumour     | KRT15             | ridge/genes    | 0.031 [-0.004, 0.062]  | -0.059 [-0.087, -0.034]       | alpha=1000.0                                    |
| tumour     | KRT17             | ridge/genes    | 0.413 [0.321, 0.510]   | -0.016 [-0.052, 0.013]        | alpha=1000.0                                    |
| tumour     | KRT6A             | ridge/genes    | 0.565 [0.247, 0.784]   | -0.008 [-0.072, 0.022]        | alpha=1000.0                                    |
| tumour     | LCN2              | ridge/genes    | 0.547 [0.461, 0.622]   | 0.007 [-0.011, 0.024]         | alpha=1000.0                                    |
| tumour     | MET               | ridge/genes    | 0.036 [0.008, 0.058]   | -0.051 [-0.073, -0.034]       | alpha=1000.0                                    |
| tumour     | OLFM4             | ridge/genes    | 0.642 [0.430, 0.759]   | 0.039 [-0.001, 0.077]         | alpha=1000.0                                    |
| tumour     | PIGR              | ridge/genes    | 0.532 [0.406, 0.629]   | -0.003 [-0.021, 0.009]        | alpha=1000.0                                    |
| tumour     | S100A2            | ridge/genes    | 0.420 [0.181, 0.604]   | -0.033 [-0.051, 0.000]        | alpha=1000.0                                    |
| tumour     | S100P             | ridge/genes    | 0.459 [0.372, 0.534]   | -0.007 [-0.049, 0.022]        | alpha=1000.0                                    |
| tumour     | SERPINB5          | ridge/genes    | 0.086 [-0.019, 0.211]  | -0.061 [-0.112, -0.027]       | alpha=1000.0                                    |
| tumour     | SLC2A1            | ridge/genes    | 0.191 [0.128, 0.249]   | -0.027 [-0.050, -0.007]       | alpha=1000.0                                    |
| tumour     | SPINK1            | ridge/genes    | 0.595 [0.484, 0.669]   | 0.016 [0.000, 0.035]          | alpha=1000.0                                    |
| tumour     | genes_pooled      | ridge/genes    | 0.412 [0.341, 0.491]   | -0.013 [-0.029, 0.002]        | alpha=1000.0                                    |
| tumour     | tumour:basal_like | ridge/genes    | 0.390 [0.287, 0.489]   | -0.008 [-0.033, 0.015]        | alpha=1000.0                                    |
| tumour     | tumour:classical  | ridge/genes    | 0.668 [0.560, 0.751]   | 0.050 [0.028, 0.071]          | alpha=1000.0                                    |
| tumour     | AGR2              | ridge/own      | 0.451 [0.362, 0.520]   | 0.000 [0.000, 0.000]          | alpha=1000.0                                    |
| tumour     | AREG              | ridge/own      | 0.157 [0.059, 0.263]   | 0.000 [0.000, 0.000]          | alpha=1000.0                                    |
| tumour     | CEACAM1           | ridge/own      | 0.132 [0.095, 0.173]   | 0.000 [0.000, 0.000]          | alpha=1000.0                                    |
| tumour     | CEACAM6           | ridge/own      | 0.627 [0.480, 0.722]   | 0.000 [0.000, 0.000]          | alpha=1000.0                                    |
| tumour     | ITGA3             | ridge/own      | 0.238 [0.201, 0.295]   | 0.000 [0.000, 0.000]          | alpha=1000.0                                    |
| tumour     | KRT14             | ridge/own      | 0.307 [0.117, 0.519]   | 0.000 [0.000, 0.000]          | alpha=1000.0                                    |
| tumour     | KRT15             | ridge/own      | 0.089 [0.072, 0.109]   | 0.000 [0.000, 0.000]          | alpha=1000.0                                    |
| tumour     | KRT17             | ridge/own      | 0.428 [0.328, 0.548]   | 0.000 [0.000, 0.000]          | alpha=1000.0                                    |
| tumour     | KRT6A             | ridge/own      | 0.574 [0.290, 0.781]   | 0.000 [0.000, 0.000]          | alpha=1000.0                                    |
| tumour     | LCN2              | ridge/own      | 0.540 [0.447, 0.622]   | 0.000 [0.000, 0.000]          | alpha=1000.0                                    |
| tumour     | MET               | ridge/own      | 0.086 [0.065, 0.106]   | 0.000 [0.000, 0.000]          | alpha=1000.0                                    |
| tumour     | OLFM4             | ridge/own      | 0.604 [0.387, 0.736]   | 0.000 [0.000, 0.000]          | alpha=1000.0                                    |
| tumour     | PIGR              | ridge/own      | 0.535 [0.411, 0.628]   | 0.000 [0.000, 0.000]          | alpha=1000.0                                    |
| tumour     | S100A2            | ridge/own      | 0.453 [0.197, 0.648]   | 0.000 [0.000, 0.000]          | alpha=1000.0                                    |
| tumour     | S100P             | ridge/own      | 0.466 [0.358, 0.568]   | 0.000 [0.000, 0.000]          | alpha=1000.0                                    |
| tumour     | SERPINB5          | ridge/own      | 0.147 [0.083, 0.255]   | 0.000 [0.000, 0.000]          | alpha=1000.0                                    |
| tumour     | SLC2A1            | ridge/own      | 0.217 [0.167, 0.269]   | 0.000 [0.000, 0.000]          | alpha=1000.0                                    |
| tumour     | SPINK1            | ridge/own      | 0.579 [0.466, 0.656]   | 0.000 [0.000, 0.000]          | alpha=1000.0                                    |
| tumour     | genes_pooled      | ridge/own      | 0.425 [0.352, 0.514]   | 0.000 [0.000, 0.000]          | alpha=1000.0                                    |
| tumour     | tumour:basal_like | ridge/own      | 0.397 [0.289, 0.512]   | 0.000 [0.000, 0.000]          | alpha=1000.0                                    |
| tumour     | tumour:classical  | ridge/own      | 0.618 [0.502, 0.711]   | 0.000 [0.000, 0.000]          | alpha=1000.0                                    |
| tumour     | AGR2              | ridge/programs | 0.461 [0.370, 0.532]   | 0.010 [0.003, 0.016]          | alpha=100.0                                     |
| tumour     | AREG              | ridge/programs | 0.149 [0.051, 0.253]   | -0.008 [-0.013, -0.003]       | alpha=100.0                                     |
| tumour     | CEACAM1           | ridge/programs | 0.138 [0.098, 0.183]   | 0.005 [-0.001, 0.011]         | alpha=100.0                                     |
| tumour     | CEACAM6           | ridge/programs | 0.654 [0.516, 0.743]   | 0.027 [0.018, 0.039]          | alpha=100.0                                     |
| tumour     | ITGA3             | ridge/programs | 0.239 [0.198, 0.297]   | 0.001 [-0.006, 0.007]         | alpha=100.0                                     |
| tumour     | KRT14             | ridge/programs | 0.301 [0.097, 0.527]   | -0.006 [-0.029, 0.012]        | alpha=100.0                                     |
| tumour     | KRT15             | ridge/programs | 0.079 [0.061, 0.099]   | -0.010 [-0.018, -0.004]       | alpha=100.0                                     |
| tumour     | KRT17             | ridge/programs | 0.433 [0.335, 0.551]   | 0.005 [-0.003, 0.014]         | alpha=100.0                                     |
| tumour     | KRT6A             | ridge/programs | 0.588 [0.312, 0.792]   | 0.015 [0.002, 0.030]          | alpha=100.0                                     |
| tumour     | LCN2              | ridge/programs | 0.556 [0.465, 0.637]   | 0.017 [0.006, 0.027]          | alpha=100.0                                     |
| tumour     | MET               | ridge/programs | 0.090 [0.070, 0.109]   | 0.004 [0.000, 0.008]          | alpha=100.0                                     |
| tumour     | OLFM4             | ridge/programs | 0.616 [0.399, 0.747]   | 0.013 [-0.004, 0.025]         | alpha=100.0                                     |
| tumour     | PIGR              | ridge/programs | 0.549 [0.426, 0.643]   | 0.014 [0.006, 0.021]          | alpha=100.0                                     |
| tumour     | S100A2            | ridge/programs | 0.456 [0.204, 0.659]   | 0.004 [-0.009, 0.016]         | alpha=100.0                                     |
| tumour     | S100P             | ridge/programs | 0.485 [0.381, 0.582]   | 0.020 [0.011, 0.027]          | alpha=100.0                                     |
| tumour     | SERPINB5          | ridge/programs | 0.153 [0.085, 0.267]   | 0.006 [0.001, 0.012]          | alpha=100.0                                     |
| tumour     | SLC2A1            | ridge/programs | 0.224 [0.170, 0.276]   | 0.007 [-0.001, 0.012]         | alpha=100.0                                     |
| tumour     | SPINK1            | ridge/programs | 0.593 [0.479, 0.670]   | 0.014 [0.005, 0.021]          | alpha=100.0                                     |
| tumour     | genes_pooled      | ridge/programs | 0.435 [0.361, 0.524]   | 0.010 [0.006, 0.014]          | alpha=100.0                                     |
| tumour     | tumour:basal_like | ridge/programs | 0.438 [0.325, 0.552]   | 0.041 [0.029, 0.052]          | alpha=100.0                                     |
| tumour     | tumour:classical  | ridge/programs | 0.693 [0.597, 0.766]   | 0.074 [0.052, 0.105]          | alpha=100.0                                     |
| tumour     | AGR2              | xgb/genes      | 0.000 [0.000, 0.000]   | -0.451 [-0.520, -0.362]       | lr=0.05,leaves=15                               |
| tumour     | AREG              | xgb/genes      | 0.000 [0.000, 0.000]   | -0.157 [-0.263, -0.059]       | lr=0.05,leaves=15                               |
| tumour     | CEACAM1           | xgb/genes      | 0.000 [0.000, 0.000]   | -0.132 [-0.173, -0.095]       | lr=0.05,leaves=15                               |
| tumour     | CEACAM6           | xgb/genes      | 0.000 [0.000, 0.000]   | -0.627 [-0.722, -0.480]       | lr=0.05,leaves=15                               |
| tumour     | ITGA3             | xgb/genes      | 0.000 [0.000, 0.000]   | -0.238 [-0.295, -0.201]       | lr=0.05,leaves=15                               |
| tumour     | KRT14             | xgb/genes      | 0.000 [0.000, 0.000]   | -0.307 [-0.519, -0.117]       | lr=0.05,leaves=15                               |
| tumour     | KRT15             | xgb/genes      | 0.000 [0.000, 0.000]   | -0.089 [-0.109, -0.072]       | lr=0.05,leaves=15                               |
| tumour     | KRT17             | xgb/genes      | 0.000 [0.000, 0.000]   | -0.428 [-0.548, -0.328]       | lr=0.05,leaves=15                               |
| tumour     | KRT6A             | xgb/genes      | 0.000 [0.000, 0.000]   | -0.574 [-0.781, -0.290]       | lr=0.05,leaves=15                               |
| tumour     | LCN2              | xgb/genes      | 0.000 [0.000, 0.000]   | -0.540 [-0.622, -0.447]       | lr=0.05,leaves=15                               |
| tumour     | MET               | xgb/genes      | 0.000 [0.000, 0.000]   | -0.086 [-0.106, -0.065]       | lr=0.05,leaves=15                               |
| tumour     | OLFM4             | xgb/genes      | 0.000 [0.000, 0.000]   | -0.604 [-0.736, -0.387]       | lr=0.05,leaves=15                               |
| tumour     | PIGR              | xgb/genes      | 0.000 [0.000, 0.000]   | -0.535 [-0.628, -0.411]       | lr=0.05,leaves=15                               |
| tumour     | S100A2            | xgb/genes      | 0.000 [0.000, 0.000]   | -0.453 [-0.648, -0.197]       | lr=0.05,leaves=15                               |
| tumour     | S100P             | xgb/genes      | 0.000 [0.000, 0.000]   | -0.466 [-0.568, -0.358]       | lr=0.05,leaves=15                               |
| tumour     | SERPINB5          | xgb/genes      | 0.000 [0.000, 0.000]   | -0.147 [-0.255, -0.083]       | lr=0.05,leaves=15                               |
| tumour     | SLC2A1            | xgb/genes      | 0.000 [0.000, 0.000]   | -0.217 [-0.269, -0.167]       | lr=0.05,leaves=15                               |
| tumour     | SPINK1            | xgb/genes      | 0.000 [0.000, 0.000]   | -0.579 [-0.656, -0.466]       | lr=0.05,leaves=15                               |
| tumour     | genes_pooled      | xgb/genes      | 0.000 [0.000, 0.000]   | -0.425 [-0.514, -0.352]       | lr=0.05,leaves=15                               |
| tumour     | tumour:basal_like | xgb/genes      | 0.444 [0.341, 0.550]   | 0.046 [0.033, 0.061]          | lr=0.05,leaves=15                               |
| tumour     | tumour:classical  | xgb/genes      | 0.694 [0.601, 0.768]   | 0.076 [0.054, 0.107]          | lr=0.05,leaves=15                               |
| tumour     | AGR2              | xgb/programs   | 0.000 [0.000, 0.000]   | -0.451 [-0.520, -0.362]       | lr=0.05,leaves=15                               |
| tumour     | AREG              | xgb/programs   | 0.000 [0.000, 0.000]   | -0.157 [-0.263, -0.059]       | lr=0.05,leaves=15                               |
| tumour     | CEACAM1           | xgb/programs   | 0.000 [0.000, 0.000]   | -0.132 [-0.173, -0.095]       | lr=0.05,leaves=15                               |
| tumour     | CEACAM6           | xgb/programs   | 0.000 [0.000, 0.000]   | -0.627 [-0.722, -0.480]       | lr=0.05,leaves=15                               |
| tumour     | ITGA3             | xgb/programs   | 0.000 [0.000, 0.000]   | -0.238 [-0.295, -0.201]       | lr=0.05,leaves=15                               |
| tumour     | KRT14             | xgb/programs   | 0.000 [0.000, 0.000]   | -0.307 [-0.519, -0.117]       | lr=0.05,leaves=15                               |
| tumour     | KRT15             | xgb/programs   | 0.000 [0.000, 0.000]   | -0.089 [-0.109, -0.072]       | lr=0.05,leaves=15                               |
| tumour     | KRT17             | xgb/programs   | 0.000 [0.000, 0.000]   | -0.428 [-0.548, -0.328]       | lr=0.05,leaves=15                               |
| tumour     | KRT6A             | xgb/programs   | 0.000 [0.000, 0.000]   | -0.574 [-0.781, -0.290]       | lr=0.05,leaves=15                               |
| tumour     | LCN2              | xgb/programs   | 0.000 [0.000, 0.000]   | -0.540 [-0.622, -0.447]       | lr=0.05,leaves=15                               |
| tumour     | MET               | xgb/programs   | 0.000 [0.000, 0.000]   | -0.086 [-0.106, -0.065]       | lr=0.05,leaves=15                               |
| tumour     | OLFM4             | xgb/programs   | 0.000 [0.000, 0.000]   | -0.604 [-0.736, -0.387]       | lr=0.05,leaves=15                               |
| tumour     | PIGR              | xgb/programs   | 0.000 [0.000, 0.000]   | -0.535 [-0.628, -0.411]       | lr=0.05,leaves=15                               |
| tumour     | S100A2            | xgb/programs   | 0.000 [0.000, 0.000]   | -0.453 [-0.648, -0.197]       | lr=0.05,leaves=15                               |
| tumour     | S100P             | xgb/programs   | 0.000 [0.000, 0.000]   | -0.466 [-0.568, -0.358]       | lr=0.05,leaves=15                               |
| tumour     | SERPINB5          | xgb/programs   | 0.000 [0.000, 0.000]   | -0.147 [-0.255, -0.083]       | lr=0.05,leaves=15                               |
| tumour     | SLC2A1            | xgb/programs   | 0.000 [0.000, 0.000]   | -0.217 [-0.269, -0.167]       | lr=0.05,leaves=15                               |
| tumour     | SPINK1            | xgb/programs   | 0.000 [0.000, 0.000]   | -0.579 [-0.656, -0.466]       | lr=0.05,leaves=15                               |
| tumour     | genes_pooled      | xgb/programs   | 0.000 [0.000, 0.000]   | -0.425 [-0.514, -0.352]       | lr=0.05,leaves=15                               |
| tumour     | tumour:basal_like | xgb/programs   | 0.438 [0.329, 0.551]   | 0.040 [0.030, 0.051]          | lr=0.05,leaves=15                               |
| tumour     | tumour:classical  | xgb/programs   | 0.693 [0.599, 0.766]   | 0.074 [0.051, 0.106]          | lr=0.05,leaves=15                               |

## Phase 3 graph models

| receiver   | target            | arch      |     lr | R2 [CI]              | per seed             | strongest baseline     | gain vs baseline [CI]   | beats   |
|:-----------|:------------------|:----------|-------:|:---------------------|:---------------------|:-----------------------|:------------------------|:--------|
| myeloid    | C1QA              | sage      | 0.0003 | 0.192 [0.117, 0.237] | 0.1921;0.1938;0.1916 | ridge/programs (0.197) | -0.004 [-0.013, 0.005]  | False   |
| myeloid    | C1QB              | sage      | 0.0003 | 0.279 [0.197, 0.340] | 0.2785;0.2766;0.2830 | ridge/programs (0.261) | 0.019 [-0.004, 0.041]   | False   |
| myeloid    | C1QC              | sage      | 0.0003 | 0.265 [0.164, 0.342] | 0.2645;0.2653;0.2637 | ridge/programs (0.251) | 0.014 [-0.006, 0.030]   | False   |
| myeloid    | CD163             | sage      | 0.0003 | 0.195 [0.108, 0.255] | 0.1953;0.1946;0.1964 | ridge/programs (0.190) | 0.006 [-0.014, 0.022]   | False   |
| myeloid    | CD9               | sage      | 0.0003 | 0.124 [0.087, 0.149] | 0.1234;0.1254;0.1234 | ridge/programs (0.112) | 0.012 [0.008, 0.018]    | True    |
| myeloid    | FABP5             | sage      | 0.0003 | 0.141 [0.095, 0.174] | 0.1423;0.1391;0.1419 | ridge/programs (0.124) | 0.017 [0.011, 0.024]    | True    |
| myeloid    | GPNMB             | sage      | 0.0003 | 0.228 [0.172, 0.269] | 0.2259;0.2293;0.2285 | ridge/programs (0.212) | 0.016 [0.009, 0.022]    | True    |
| myeloid    | LGALS3            | sage      | 0.0003 | 0.152 [0.124, 0.178] | 0.1501;0.1522;0.1550 | ridge/programs (0.142) | 0.010 [0.007, 0.015]    | True    |
| myeloid    | MRC1              | sage      | 0.0003 | 0.092 [0.070, 0.107] | 0.0906;0.0943;0.0916 | ridge/programs (0.095) | -0.003 [-0.014, 0.008]  | False   |
| myeloid    | SELENOP           | sage      | 0.0003 | 0.271 [0.208, 0.317] | 0.2704;0.2680;0.2755 | ridge/programs (0.245) | 0.026 [0.013, 0.044]    | True    |
| myeloid    | SPP1              | sage      | 0.0003 | 0.654 [0.519, 0.745] | 0.6535;0.6538;0.6541 | ridge/programs (0.633) | 0.020 [0.012, 0.038]    | True    |
| myeloid    | genes_pooled      | sage      | 0.0003 | 0.289 [0.214, 0.348] | 0.2888;0.2893;0.2903 | ridge/programs (0.277) | 0.013 [0.005, 0.020]    | True    |
| myeloid    | myeloid:C1Q_TAM   | sage      | 0.0003 | 0.336 [0.220, 0.432] | 0.3296;0.3401;0.3386 | xgb/genes (0.334)      | 0.003 [-0.013, 0.015]   | False   |
| myeloid    | myeloid:SPP1_TAM  | sage      | 0.0003 | 0.377 [0.301, 0.428] | 0.3754;0.3735;0.3828 | xgb/programs (0.387)   | -0.010 [-0.021, 0.007]  | False   |
| myeloid    | C1QA              | attention | 0.0003 | 0.191 [0.115, 0.236] | 0.1916;0.1883;0.1938 | ridge/programs (0.197) | -0.005 [-0.011, 0.003]  | False   |
| myeloid    | C1QB              | attention | 0.0003 | 0.273 [0.187, 0.336] | 0.2743;0.2712;0.2735 | ridge/programs (0.261) | 0.012 [-0.012, 0.037]   | False   |
| myeloid    | C1QC              | attention | 0.0003 | 0.264 [0.161, 0.342] | 0.2647;0.2612;0.2654 | ridge/programs (0.251) | 0.013 [-0.006, 0.030]   | False   |
| myeloid    | CD163             | attention | 0.0003 | 0.191 [0.104, 0.251] | 0.1909;0.1884;0.1936 | ridge/programs (0.190) | 0.001 [-0.020, 0.020]   | False   |
| myeloid    | CD9               | attention | 0.0003 | 0.120 [0.082, 0.146] | 0.1251;0.1186;0.1174 | ridge/programs (0.112) | 0.008 [0.004, 0.012]    | True    |
| myeloid    | FABP5             | attention | 0.0003 | 0.137 [0.093, 0.168] | 0.1371;0.1371;0.1365 | ridge/programs (0.124) | 0.013 [0.009, 0.019]    | True    |
| myeloid    | GPNMB             | attention | 0.0003 | 0.222 [0.171, 0.264] | 0.2234;0.2196;0.2242 | ridge/programs (0.212) | 0.010 [0.002, 0.020]    | True    |
| myeloid    | LGALS3            | attention | 0.0003 | 0.152 [0.121, 0.180] | 0.1452;0.1572;0.1537 | ridge/programs (0.142) | 0.010 [0.006, 0.014]    | True    |
| myeloid    | MRC1              | attention | 0.0003 | 0.085 [0.059, 0.103] | 0.0834;0.0833;0.0878 | ridge/programs (0.095) | -0.010 [-0.018, -0.002] | False   |
| myeloid    | SELENOP           | attention | 0.0003 | 0.260 [0.197, 0.308] | 0.2634;0.2535;0.2638 | ridge/programs (0.245) | 0.015 [-0.004, 0.036]   | False   |
| myeloid    | SPP1              | attention | 0.0003 | 0.656 [0.520, 0.746] | 0.6556;0.6555;0.6556 | ridge/programs (0.633) | 0.022 [0.013, 0.040]    | True    |
| myeloid    | genes_pooled      | attention | 0.0003 | 0.286 [0.210, 0.344] | 0.2865;0.2847;0.2875 | ridge/programs (0.277) | 0.010 [0.002, 0.018]    | True    |
| myeloid    | myeloid:C1Q_TAM   | attention | 0.0003 | 0.335 [0.215, 0.431] | 0.3340;0.3320;0.3377 | xgb/genes (0.334)      | 0.001 [-0.016, 0.014]   | False   |
| myeloid    | myeloid:SPP1_TAM  | attention | 0.0003 | 0.381 [0.307, 0.432] | 0.3798;0.3776;0.3858 | xgb/programs (0.387)   | -0.006 [-0.017, 0.008]  | False   |
| tumour     | AGR2              | sage      | 0.0003 | 0.465 [0.366, 0.544] | 0.4588;0.4694;0.4664 | ridge/programs (0.461) | 0.004 [-0.016, 0.019]   | False   |
| tumour     | AREG              | sage      | 0.0003 | 0.146 [0.060, 0.234] | 0.1521;0.1371;0.1499 | ridge/own (0.157)      | -0.011 [-0.028, 0.006]  | False   |
| tumour     | CEACAM1           | sage      | 0.0003 | 0.139 [0.092, 0.194] | 0.1337;0.1414;0.1417 | ridge/programs (0.138) | 0.001 [-0.012, 0.012]   | False   |
| tumour     | CEACAM6           | sage      | 0.0003 | 0.651 [0.522, 0.743] | 0.6465;0.6569;0.6504 | ridge/programs (0.654) | -0.002 [-0.032, 0.018]  | False   |
| tumour     | ITGA3             | sage      | 0.0003 | 0.227 [0.189, 0.274] | 0.2314;0.2211;0.2279 | ridge/programs (0.239) | -0.012 [-0.032, 0.008]  | False   |
| tumour     | KRT14             | sage      | 0.0003 | 0.307 [0.116, 0.523] | 0.3102;0.3077;0.3017 | ridge/own (0.307)      | -0.000 [-0.012, 0.010]  | False   |
| tumour     | KRT15             | sage      | 0.0003 | 0.083 [0.062, 0.107] | 0.0876;0.0789;0.0834 | ridge/own (0.089)      | -0.006 [-0.024, 0.007]  | False   |
| tumour     | KRT17             | sage      | 0.0003 | 0.413 [0.336, 0.504] | 0.4199;0.3982;0.4212 | ridge/programs (0.433) | -0.020 [-0.049, 0.009]  | False   |
| tumour     | KRT6A             | sage      | 0.0003 | 0.590 [0.297, 0.799] | 0.5860;0.5940;0.5901 | ridge/programs (0.588) | 0.002 [-0.025, 0.019]   | False   |
| tumour     | LCN2              | sage      | 0.0003 | 0.554 [0.465, 0.632] | 0.5510;0.5583;0.5526 | ridge/programs (0.556) | -0.003 [-0.013, 0.006]  | False   |
| tumour     | MET               | sage      | 0.0003 | 0.076 [0.053, 0.098] | 0.0764;0.0731;0.0788 | ridge/programs (0.090) | -0.014 [-0.028, -0.003] | False   |
| tumour     | OLFM4             | sage      | 0.0003 | 0.648 [0.443, 0.762] | 0.6472;0.6500;0.6466 | ridge/genes (0.642)    | 0.006 [-0.008, 0.024]   | False   |
| tumour     | PIGR              | sage      | 0.0003 | 0.554 [0.431, 0.647] | 0.5546;0.5530;0.5551 | ridge/programs (0.549) | 0.005 [-0.000, 0.010]   | False   |
| tumour     | S100A2            | sage      | 0.0003 | 0.442 [0.207, 0.675] | 0.4368;0.4393;0.4495 | ridge/programs (0.456) | -0.014 [-0.084, 0.028]  | False   |
| tumour     | S100P             | sage      | 0.0003 | 0.474 [0.377, 0.566] | 0.4774;0.4715;0.4737 | ridge/programs (0.485) | -0.011 [-0.027, 0.005]  | False   |
| tumour     | SERPINB5          | sage      | 0.0003 | 0.148 [0.075, 0.260] | 0.1491;0.1459;0.1501 | ridge/programs (0.153) | -0.005 [-0.022, 0.009]  | False   |
| tumour     | SLC2A1            | sage      | 0.0003 | 0.236 [0.178, 0.293] | 0.2351;0.2363;0.2377 | ridge/programs (0.224) | 0.012 [0.002, 0.022]    | True    |
| tumour     | SPINK1            | sage      | 0.0003 | 0.603 [0.489, 0.679] | 0.6047;0.6037;0.5992 | ridge/genes (0.595)    | 0.007 [0.001, 0.015]    | True    |
| tumour     | genes_pooled      | sage      | 0.0003 | 0.436 [0.364, 0.523] | 0.4356;0.4355;0.4364 | ridge/programs (0.435) | 0.000 [-0.005, 0.007]   | False   |
| tumour     | tumour:basal_like | sage      | 0.0003 | 0.281 [0.196, 0.367] | 0.2661;0.2940;0.2833 | xgb/genes (0.444)      | -0.163 [-0.330, -0.045] | False   |
| tumour     | tumour:classical  | sage      | 0.0003 | 0.597 [0.444, 0.703] | 0.5959;0.5973;0.5991 | xgb/genes (0.694)      | -0.097 [-0.172, -0.053] | False   |
| tumour     | AGR2              | attention | 0.0003 | 0.450 [0.346, 0.533] | 0.4508;0.4489;0.4491 | ridge/programs (0.461) | -0.011 [-0.037, 0.011]  | False   |
| tumour     | AREG              | attention | 0.0003 | 0.134 [0.055, 0.207] | 0.1348;0.1430;0.1237 | ridge/own (0.157)      | -0.024 [-0.056, 0.001]  | False   |
| tumour     | CEACAM1           | attention | 0.0003 | 0.134 [0.089, 0.184] | 0.1366;0.1295;0.1364 | ridge/programs (0.138) | -0.004 [-0.016, 0.005]  | False   |
| tumour     | CEACAM6           | attention | 0.0003 | 0.643 [0.512, 0.734] | 0.6312;0.6549;0.6421 | ridge/programs (0.654) | -0.011 [-0.037, 0.009]  | False   |
| tumour     | ITGA3             | attention | 0.0003 | 0.218 [0.178, 0.268] | 0.2188;0.2122;0.2239 | ridge/programs (0.239) | -0.021 [-0.036, -0.005] | False   |
| tumour     | KRT14             | attention | 0.0003 | 0.305 [0.116, 0.518] | 0.3082;0.2979;0.3100 | ridge/own (0.307)      | -0.001 [-0.012, 0.007]  | False   |
| tumour     | KRT15             | attention | 0.0003 | 0.081 [0.063, 0.102] | 0.0865;0.0804;0.0769 | ridge/own (0.089)      | -0.008 [-0.025, 0.003]  | False   |
| tumour     | KRT17             | attention | 0.0003 | 0.414 [0.320, 0.514] | 0.4251;0.4003;0.4155 | ridge/programs (0.433) | -0.020 [-0.043, 0.002]  | False   |
| tumour     | KRT6A             | attention | 0.0003 | 0.596 [0.298, 0.806] | 0.5924;0.5978;0.5987 | ridge/programs (0.588) | 0.008 [-0.021, 0.025]   | False   |
| tumour     | LCN2              | attention | 0.0003 | 0.548 [0.454, 0.627] | 0.5391;0.5508;0.5527 | ridge/programs (0.556) | -0.009 [-0.018, -0.001] | False   |
| tumour     | MET               | attention | 0.0003 | 0.076 [0.055, 0.097] | 0.0781;0.0795;0.0714 | ridge/programs (0.090) | -0.014 [-0.022, -0.006] | False   |
| tumour     | OLFM4             | attention | 0.0003 | 0.656 [0.447, 0.773] | 0.6563;0.6527;0.6591 | ridge/genes (0.642)    | 0.014 [0.005, 0.026]    | True    |
| tumour     | PIGR              | attention | 0.0003 | 0.555 [0.430, 0.650] | 0.5532;0.5550;0.5564 | ridge/programs (0.549) | 0.006 [0.003, 0.008]    | True    |
| tumour     | S100A2            | attention | 0.0003 | 0.422 [0.179, 0.665] | 0.4192;0.4163;0.4297 | ridge/programs (0.456) | -0.035 [-0.124, 0.020]  | False   |
| tumour     | S100P             | attention | 0.0003 | 0.467 [0.367, 0.559] | 0.4676;0.4681;0.4666 | ridge/programs (0.485) | -0.018 [-0.030, -0.007] | False   |
| tumour     | SERPINB5          | attention | 0.0003 | 0.145 [0.073, 0.257] | 0.1498;0.1447;0.1405 | ridge/programs (0.153) | -0.008 [-0.024, 0.006]  | False   |
| tumour     | SLC2A1            | attention | 0.0003 | 0.233 [0.181, 0.287] | 0.2315;0.2318;0.2343 | ridge/programs (0.224) | 0.009 [-0.000, 0.018]   | False   |
| tumour     | SPINK1            | attention | 0.0003 | 0.598 [0.481, 0.675] | 0.5968;0.5990;0.5977 | ridge/genes (0.595)    | 0.003 [-0.007, 0.009]   | False   |
| tumour     | genes_pooled      | attention | 0.0003 | 0.432 [0.359, 0.519] | 0.4309;0.4313;0.4329 | ridge/programs (0.435) | -0.004 [-0.010, 0.003]  | False   |
| tumour     | tumour:basal_like | attention | 0.0003 | 0.295 [0.234, 0.367] | 0.2990;0.2908;0.2966 | xgb/genes (0.444)      | -0.148 [-0.279, -0.051] | False   |
| tumour     | tumour:classical  | attention | 0.0003 | 0.601 [0.448, 0.710] | 0.5851;0.5940;0.6237 | xgb/genes (0.694)      | -0.093 [-0.165, -0.051] | False   |

## Phase 3 ablations (selected architecture)

| receiver   | target            | ablation    | R2 ablation [CI]     | full - ablation [CI]    | CI excludes 0   | quantity                                                |
|:-----------|:------------------|:------------|:---------------------|:------------------------|:----------------|:--------------------------------------------------------|
| myeloid    | C1QA              | same_type   | 0.183 [0.107, 0.228] | 0.010 [0.001, 0.019]    | True            | cross-type contribution (full - same-type only)         |
| myeloid    | C1QA              | other_type  | 0.175 [0.107, 0.215] | 0.018 [0.007, 0.024]    | True            | full - other-type only                                  |
| myeloid    | C1QA              | composition | 0.178 [0.105, 0.221] | 0.015 [0.004, 0.025]    | True            | full - composition only                                 |
| myeloid    | C1QA              | shuffled    | 0.176 [0.103, 0.220] | 0.017 [0.011, 0.022]    | True            | full - shuffled neighbour expression (zonation control) |
| myeloid    | C1QA              | no_graph    | 0.180 [0.103, 0.226] | 0.013 [0.001, 0.025]    | True            | full - no graph                                         |
| myeloid    | C1QB              | same_type   | 0.271 [0.186, 0.329] | 0.009 [0.001, 0.016]    | True            | cross-type contribution (full - same-type only)         |
| myeloid    | C1QB              | other_type  | 0.258 [0.183, 0.312] | 0.021 [0.007, 0.030]    | True            | full - other-type only                                  |
| myeloid    | C1QB              | composition | 0.262 [0.175, 0.318] | 0.017 [0.001, 0.034]    | True            | full - composition only                                 |
| myeloid    | C1QB              | shuffled    | 0.255 [0.173, 0.313] | 0.025 [0.015, 0.032]    | True            | full - shuffled neighbour expression (zonation control) |
| myeloid    | C1QB              | no_graph    | 0.261 [0.175, 0.317] | 0.018 [0.000, 0.036]    | True            | full - no graph                                         |
| myeloid    | C1QC              | same_type   | 0.253 [0.151, 0.327] | 0.011 [-0.001, 0.021]   | False           | cross-type contribution (full - same-type only)         |
| myeloid    | C1QC              | other_type  | 0.244 [0.146, 0.318] | 0.020 [0.011, 0.027]    | True            | full - other-type only                                  |
| myeloid    | C1QC              | composition | 0.245 [0.144, 0.311] | 0.019 [-0.002, 0.039]   | False           | full - composition only                                 |
| myeloid    | C1QC              | shuffled    | 0.237 [0.139, 0.309] | 0.027 [0.015, 0.037]    | True            | full - shuffled neighbour expression (zonation control) |
| myeloid    | C1QC              | no_graph    | 0.243 [0.138, 0.313] | 0.021 [0.001, 0.040]    | True            | full - no graph                                         |
| myeloid    | CD163             | same_type   | 0.190 [0.104, 0.247] | 0.006 [-0.005, 0.014]   | False           | cross-type contribution (full - same-type only)         |
| myeloid    | CD163             | other_type  | 0.183 [0.096, 0.241] | 0.013 [0.002, 0.019]    | True            | full - other-type only                                  |
| myeloid    | CD163             | composition | 0.183 [0.104, 0.239] | 0.012 [-0.005, 0.024]   | False           | full - composition only                                 |
| myeloid    | CD163             | shuffled    | 0.184 [0.094, 0.243] | 0.012 [0.007, 0.016]    | True            | full - shuffled neighbour expression (zonation control) |
| myeloid    | CD163             | no_graph    | 0.184 [0.107, 0.236] | 0.012 [-0.008, 0.026]   | False           | full - no graph                                         |
| myeloid    | CD9               | same_type   | 0.117 [0.081, 0.142] | 0.007 [0.004, 0.009]    | True            | cross-type contribution (full - same-type only)         |
| myeloid    | CD9               | other_type  | 0.109 [0.071, 0.136] | 0.015 [0.011, 0.019]    | True            | full - other-type only                                  |
| myeloid    | CD9               | composition | 0.114 [0.074, 0.140] | 0.010 [0.006, 0.017]    | True            | full - composition only                                 |
| myeloid    | CD9               | shuffled    | 0.113 [0.078, 0.138] | 0.011 [0.007, 0.013]    | True            | full - shuffled neighbour expression (zonation control) |
| myeloid    | CD9               | no_graph    | 0.113 [0.076, 0.138] | 0.011 [0.008, 0.016]    | True            | full - no graph                                         |
| myeloid    | FABP5             | same_type   | 0.131 [0.091, 0.161] | 0.010 [0.003, 0.016]    | True            | cross-type contribution (full - same-type only)         |
| myeloid    | FABP5             | other_type  | 0.127 [0.087, 0.155] | 0.014 [0.008, 0.019]    | True            | full - other-type only                                  |
| myeloid    | FABP5             | composition | 0.127 [0.089, 0.155] | 0.014 [0.006, 0.020]    | True            | full - composition only                                 |
| myeloid    | FABP5             | shuffled    | 0.128 [0.085, 0.158] | 0.013 [0.008, 0.016]    | True            | full - shuffled neighbour expression (zonation control) |
| myeloid    | FABP5             | no_graph    | 0.121 [0.084, 0.147] | 0.020 [0.010, 0.028]    | True            | full - no graph                                         |
| myeloid    | GPNMB             | same_type   | 0.210 [0.162, 0.245] | 0.018 [0.009, 0.026]    | True            | cross-type contribution (full - same-type only)         |
| myeloid    | GPNMB             | other_type  | 0.193 [0.147, 0.227] | 0.035 [0.021, 0.044]    | True            | full - other-type only                                  |
| myeloid    | GPNMB             | composition | 0.206 [0.158, 0.241] | 0.022 [0.011, 0.032]    | True            | full - composition only                                 |
| myeloid    | GPNMB             | shuffled    | 0.197 [0.150, 0.234] | 0.031 [0.022, 0.038]    | True            | full - shuffled neighbour expression (zonation control) |
| myeloid    | GPNMB             | no_graph    | 0.197 [0.151, 0.229] | 0.031 [0.019, 0.042]    | True            | full - no graph                                         |
| myeloid    | LGALS3            | same_type   | 0.147 [0.119, 0.172] | 0.006 [0.002, 0.010]    | True            | cross-type contribution (full - same-type only)         |
| myeloid    | LGALS3            | other_type  | 0.138 [0.111, 0.162] | 0.015 [0.009, 0.019]    | True            | full - other-type only                                  |
| myeloid    | LGALS3            | composition | 0.145 [0.119, 0.170] | 0.007 [0.003, 0.012]    | True            | full - composition only                                 |
| myeloid    | LGALS3            | shuffled    | 0.142 [0.116, 0.167] | 0.010 [0.007, 0.013]    | True            | full - shuffled neighbour expression (zonation control) |
| myeloid    | LGALS3            | no_graph    | 0.138 [0.108, 0.165] | 0.014 [0.010, 0.017]    | True            | full - no graph                                         |
| myeloid    | MRC1              | same_type   | 0.091 [0.063, 0.111] | 0.001 [-0.006, 0.012]   | False           | cross-type contribution (full - same-type only)         |
| myeloid    | MRC1              | other_type  | 0.086 [0.062, 0.103] | 0.007 [0.002, 0.012]    | True            | full - other-type only                                  |
| myeloid    | MRC1              | composition | 0.092 [0.065, 0.109] | 0.000 [-0.004, 0.006]   | False           | full - composition only                                 |
| myeloid    | MRC1              | shuffled    | 0.090 [0.066, 0.107] | 0.002 [-0.001, 0.005]   | False           | full - shuffled neighbour expression (zonation control) |
| myeloid    | MRC1              | no_graph    | 0.088 [0.062, 0.105] | 0.004 [-0.000, 0.011]   | False           | full - no graph                                         |
| myeloid    | SELENOP           | same_type   | 0.240 [0.173, 0.285] | 0.032 [0.019, 0.052]    | True            | cross-type contribution (full - same-type only)         |
| myeloid    | SELENOP           | other_type  | 0.248 [0.189, 0.291] | 0.024 [0.015, 0.029]    | True            | full - other-type only                                  |
| myeloid    | SELENOP           | composition | 0.241 [0.176, 0.285] | 0.030 [0.016, 0.051]    | True            | full - composition only                                 |
| myeloid    | SELENOP           | shuffled    | 0.242 [0.182, 0.287] | 0.029 [0.019, 0.036]    | True            | full - shuffled neighbour expression (zonation control) |
| myeloid    | SELENOP           | no_graph    | 0.234 [0.169, 0.276] | 0.038 [0.019, 0.061]    | True            | full - no graph                                         |
| myeloid    | SPP1              | same_type   | 0.628 [0.493, 0.718] | 0.025 [0.016, 0.036]    | True            | cross-type contribution (full - same-type only)         |
| myeloid    | SPP1              | other_type  | 0.622 [0.498, 0.710] | 0.032 [0.019, 0.039]    | True            | full - other-type only                                  |
| myeloid    | SPP1              | composition | 0.621 [0.488, 0.713] | 0.033 [0.020, 0.052]    | True            | full - composition only                                 |
| myeloid    | SPP1              | shuffled    | 0.621 [0.494, 0.713] | 0.033 [0.018, 0.050]    | True            | full - shuffled neighbour expression (zonation control) |
| myeloid    | SPP1              | no_graph    | 0.615 [0.479, 0.709] | 0.039 [0.023, 0.066]    | True            | full - no graph                                         |
| myeloid    | genes_pooled      | same_type   | 0.276 [0.200, 0.332] | 0.014 [0.009, 0.018]    | True            | cross-type contribution (full - same-type only)         |
| myeloid    | genes_pooled      | other_type  | 0.268 [0.199, 0.324] | 0.021 [0.013, 0.026]    | True            | full - other-type only                                  |
| myeloid    | genes_pooled      | composition | 0.271 [0.198, 0.326] | 0.019 [0.010, 0.028]    | True            | full - composition only                                 |
| myeloid    | genes_pooled      | shuffled    | 0.268 [0.196, 0.324] | 0.021 [0.015, 0.027]    | True            | full - shuffled neighbour expression (zonation control) |
| myeloid    | genes_pooled      | no_graph    | 0.267 [0.194, 0.324] | 0.022 [0.013, 0.033]    | True            | full - no graph                                         |
| myeloid    | myeloid:C1Q_TAM   | same_type   | 0.317 [0.205, 0.410] | 0.019 [0.004, 0.028]    | True            | cross-type contribution (full - same-type only)         |
| myeloid    | myeloid:C1Q_TAM   | other_type  | 0.308 [0.199, 0.402] | 0.028 [0.013, 0.035]    | True            | full - other-type only                                  |
| myeloid    | myeloid:C1Q_TAM   | composition | 0.308 [0.194, 0.390] | 0.028 [0.004, 0.048]    | True            | full - composition only                                 |
| myeloid    | myeloid:C1Q_TAM   | shuffled    | 0.297 [0.182, 0.388] | 0.039 [0.026, 0.049]    | True            | full - shuffled neighbour expression (zonation control) |
| myeloid    | myeloid:C1Q_TAM   | no_graph    | 0.307 [0.191, 0.394] | 0.029 [0.005, 0.051]    | True            | full - no graph                                         |
| myeloid    | myeloid:SPP1_TAM  | same_type   | 0.362 [0.282, 0.422] | 0.015 [0.004, 0.028]    | True            | cross-type contribution (full - same-type only)         |
| myeloid    | myeloid:SPP1_TAM  | other_type  | 0.330 [0.254, 0.384] | 0.047 [0.038, 0.055]    | True            | full - other-type only                                  |
| myeloid    | myeloid:SPP1_TAM  | composition | 0.347 [0.268, 0.411] | 0.031 [0.014, 0.050]    | True            | full - composition only                                 |
| myeloid    | myeloid:SPP1_TAM  | shuffled    | 0.342 [0.267, 0.398] | 0.035 [0.025, 0.048]    | True            | full - shuffled neighbour expression (zonation control) |
| myeloid    | myeloid:SPP1_TAM  | no_graph    | 0.325 [0.254, 0.381] | 0.052 [0.037, 0.071]    | True            | full - no graph                                         |
| tumour     | AGR2              | same_type   | 0.469 [0.371, 0.547] | -0.004 [-0.010, 0.000]  | False           | cross-type contribution (full - same-type only)         |
| tumour     | AGR2              | other_type  | 0.431 [0.334, 0.509] | 0.034 [0.023, 0.045]    | True            | full - other-type only                                  |
| tumour     | AGR2              | no_graph    | 0.449 [0.357, 0.524] | 0.016 [-0.003, 0.033]   | False           | full - no graph                                         |
| tumour     | AREG              | same_type   | 0.140 [0.060, 0.223] | 0.006 [-0.001, 0.012]   | False           | cross-type contribution (full - same-type only)         |
| tumour     | AREG              | other_type  | 0.137 [0.056, 0.222] | 0.009 [0.003, 0.014]    | True            | full - other-type only                                  |
| tumour     | AREG              | no_graph    | 0.137 [0.059, 0.216] | 0.010 [-0.000, 0.018]   | False           | full - no graph                                         |
| tumour     | CEACAM1           | same_type   | 0.142 [0.097, 0.193] | -0.003 [-0.007, 0.001]  | False           | cross-type contribution (full - same-type only)         |
| tumour     | CEACAM1           | other_type  | 0.126 [0.089, 0.166] | 0.013 [-0.002, 0.029]   | False           | full - other-type only                                  |
| tumour     | CEACAM1           | no_graph    | 0.123 [0.085, 0.165] | 0.016 [-0.001, 0.031]   | False           | full - no graph                                         |
| tumour     | CEACAM6           | same_type   | 0.648 [0.519, 0.739] | 0.003 [-0.006, 0.010]   | False           | cross-type contribution (full - same-type only)         |
| tumour     | CEACAM6           | other_type  | 0.593 [0.444, 0.690] | 0.058 [0.033, 0.086]    | True            | full - other-type only                                  |
| tumour     | CEACAM6           | no_graph    | 0.616 [0.469, 0.712] | 0.035 [0.001, 0.065]    | True            | full - no graph                                         |
| tumour     | ITGA3             | same_type   | 0.228 [0.186, 0.283] | -0.001 [-0.011, 0.009]  | False           | cross-type contribution (full - same-type only)         |
| tumour     | ITGA3             | other_type  | 0.208 [0.161, 0.268] | 0.019 [-0.001, 0.039]   | False           | full - other-type only                                  |
| tumour     | ITGA3             | no_graph    | 0.214 [0.166, 0.274] | 0.012 [-0.009, 0.035]   | False           | full - no graph                                         |
| tumour     | KRT14             | same_type   | 0.296 [0.117, 0.499] | 0.010 [-0.003, 0.025]   | False           | cross-type contribution (full - same-type only)         |
| tumour     | KRT14             | other_type  | 0.278 [0.101, 0.474] | 0.029 [0.005, 0.051]    | True            | full - other-type only                                  |
| tumour     | KRT14             | no_graph    | 0.288 [0.107, 0.492] | 0.018 [-0.001, 0.034]   | False           | full - no graph                                         |
| tumour     | KRT15             | same_type   | 0.085 [0.069, 0.105] | -0.002 [-0.015, 0.006]  | False           | cross-type contribution (full - same-type only)         |
| tumour     | KRT15             | other_type  | 0.075 [0.060, 0.091] | 0.009 [-0.002, 0.020]   | False           | full - other-type only                                  |
| tumour     | KRT15             | no_graph    | 0.078 [0.062, 0.094] | 0.006 [-0.009, 0.018]   | False           | full - no graph                                         |
| tumour     | KRT17             | same_type   | 0.422 [0.337, 0.528] | -0.009 [-0.027, 0.007]  | False           | cross-type contribution (full - same-type only)         |
| tumour     | KRT17             | other_type  | 0.392 [0.311, 0.496] | 0.021 [-0.007, 0.050]   | False           | full - other-type only                                  |
| tumour     | KRT17             | no_graph    | 0.411 [0.325, 0.521] | 0.002 [-0.025, 0.026]   | False           | full - no graph                                         |
| tumour     | KRT6A             | same_type   | 0.597 [0.308, 0.804] | -0.007 [-0.015, -0.003] | True            | cross-type contribution (full - same-type only)         |
| tumour     | KRT6A             | other_type  | 0.512 [0.240, 0.711] | 0.078 [0.038, 0.099]    | True            | full - other-type only                                  |
| tumour     | KRT6A             | no_graph    | 0.544 [0.264, 0.751] | 0.046 [0.012, 0.071]    | True            | full - no graph                                         |
| tumour     | LCN2              | same_type   | 0.553 [0.462, 0.631] | 0.001 [-0.007, 0.008]   | False           | cross-type contribution (full - same-type only)         |
| tumour     | LCN2              | other_type  | 0.502 [0.411, 0.581] | 0.052 [0.032, 0.070]    | True            | full - other-type only                                  |
| tumour     | LCN2              | no_graph    | 0.522 [0.424, 0.605] | 0.032 [0.013, 0.053]    | True            | full - no graph                                         |
| tumour     | MET               | same_type   | 0.077 [0.057, 0.096] | -0.001 [-0.007, 0.004]  | False           | cross-type contribution (full - same-type only)         |
| tumour     | MET               | other_type  | 0.076 [0.056, 0.093] | 0.001 [-0.008, 0.008]   | False           | full - other-type only                                  |
| tumour     | MET               | no_graph    | 0.083 [0.062, 0.103] | -0.007 [-0.018, 0.001]  | False           | full - no graph                                         |
| tumour     | OLFM4             | same_type   | 0.657 [0.449, 0.773] | -0.009 [-0.016, -0.002] | True            | cross-type contribution (full - same-type only)         |
| tumour     | OLFM4             | other_type  | 0.587 [0.359, 0.722] | 0.061 [0.030, 0.097]    | True            | full - other-type only                                  |
| tumour     | OLFM4             | no_graph    | 0.595 [0.368, 0.727] | 0.053 [0.020, 0.089]    | True            | full - no graph                                         |
| tumour     | PIGR              | same_type   | 0.551 [0.429, 0.646] | 0.003 [-0.003, 0.008]   | False           | cross-type contribution (full - same-type only)         |
| tumour     | PIGR              | other_type  | 0.520 [0.395, 0.617] | 0.034 [0.023, 0.048]    | True            | full - other-type only                                  |
| tumour     | PIGR              | no_graph    | 0.529 [0.408, 0.622] | 0.025 [0.013, 0.036]    | True            | full - no graph                                         |
| tumour     | S100A2            | same_type   | 0.435 [0.199, 0.672] | 0.007 [0.002, 0.015]    | True            | cross-type contribution (full - same-type only)         |
| tumour     | S100A2            | other_type  | 0.366 [0.154, 0.562] | 0.076 [0.028, 0.114]    | True            | full - other-type only                                  |
| tumour     | S100A2            | no_graph    | 0.395 [0.179, 0.583] | 0.047 [-0.006, 0.092]   | False           | full - no graph                                         |
| tumour     | S100P             | same_type   | 0.474 [0.376, 0.566] | 0.000 [-0.006, 0.007]   | False           | cross-type contribution (full - same-type only)         |
| tumour     | S100P             | other_type  | 0.444 [0.339, 0.542] | 0.031 [0.014, 0.047]    | True            | full - other-type only                                  |
| tumour     | S100P             | no_graph    | 0.455 [0.345, 0.558] | 0.019 [-0.001, 0.042]   | False           | full - no graph                                         |
| tumour     | SERPINB5          | same_type   | 0.148 [0.076, 0.260] | 0.000 [-0.004, 0.004]   | False           | cross-type contribution (full - same-type only)         |
| tumour     | SERPINB5          | other_type  | 0.135 [0.069, 0.239] | 0.013 [0.000, 0.026]    | True            | full - other-type only                                  |
| tumour     | SERPINB5          | no_graph    | 0.143 [0.074, 0.251] | 0.006 [-0.005, 0.014]   | False           | full - no graph                                         |
| tumour     | SLC2A1            | same_type   | 0.234 [0.180, 0.288] | 0.003 [-0.006, 0.010]   | False           | cross-type contribution (full - same-type only)         |
| tumour     | SLC2A1            | other_type  | 0.209 [0.164, 0.258] | 0.027 [0.011, 0.043]    | True            | full - other-type only                                  |
| tumour     | SLC2A1            | no_graph    | 0.207 [0.166, 0.251] | 0.030 [0.008, 0.050]    | True            | full - no graph                                         |
| tumour     | SPINK1            | same_type   | 0.606 [0.493, 0.683] | -0.003 [-0.006, -0.001] | True            | cross-type contribution (full - same-type only)         |
| tumour     | SPINK1            | other_type  | 0.560 [0.433, 0.646] | 0.042 [0.029, 0.063]    | True            | full - other-type only                                  |
| tumour     | SPINK1            | no_graph    | 0.577 [0.457, 0.660] | 0.025 [0.014, 0.043]    | True            | full - no graph                                         |
| tumour     | genes_pooled      | same_type   | 0.437 [0.364, 0.524] | -0.001 [-0.002, 0.001]  | False           | cross-type contribution (full - same-type only)         |
| tumour     | genes_pooled      | other_type  | 0.397 [0.326, 0.478] | 0.039 [0.029, 0.050]    | True            | full - other-type only                                  |
| tumour     | genes_pooled      | no_graph    | 0.411 [0.338, 0.495] | 0.025 [0.016, 0.034]    | True            | full - no graph                                         |
| tumour     | tumour:basal_like | same_type   | 0.322 [0.266, 0.374] | -0.041 [-0.121, 0.011]  | False           | cross-type contribution (full - same-type only)         |
| tumour     | tumour:basal_like | other_type  | 0.319 [0.223, 0.414] | -0.038 [-0.193, 0.066]  | False           | full - other-type only                                  |
| tumour     | tumour:basal_like | no_graph    | 0.340 [0.250, 0.433] | -0.059 [-0.215, 0.045]  | False           | full - no graph                                         |
| tumour     | tumour:classical  | same_type   | 0.591 [0.424, 0.705] | 0.007 [-0.004, 0.024]   | False           | cross-type contribution (full - same-type only)         |
| tumour     | tumour:classical  | other_type  | 0.549 [0.389, 0.665] | 0.048 [0.027, 0.069]    | True            | full - other-type only                                  |
| tumour     | tumour:classical  | no_graph    | 0.567 [0.414, 0.677] | 0.031 [0.011, 0.052]    | True            | full - no graph                                         |

## Phase 4 external test (GSE310352; single touch; units UNVERIFIED as independent patients)

| receiver   | target            | model          |   internal R2 | external R2 [unit CI]   |   drop | gain vs own (ext) [CI]   | per slide                                                                                                                                                  |
|:-----------|:------------------|:---------------|--------------:|:------------------------|-------:|:-------------------------|:-----------------------------------------------------------------------------------------------------------------------------------------------------------|
| myeloid    | myeloid:C1Q_TAM   | mean           |         0     | 0.000 [0.000, 0.000]    |  0     | -0.310 [-0.463, -0.211]  | {"slide1": 0.0, "slide2": 0.0, "slide3": 0.0, "slide4": 0.0, "slide5": 0.0, "slide5b2": 0.0, "slide5b3": 0.0, "slide6": 0.0}                               |
| myeloid    | myeloid:SPP1_TAM  | mean           |         0     | 0.000 [0.000, 0.000]    |  0     | -0.325 [-0.438, -0.203]  | {"slide1": 0.0, "slide2": 0.0, "slide3": 0.0, "slide4": 0.0, "slide5": 0.0, "slide5b2": 0.0, "slide5b3": 0.0, "slide6": 0.0}                               |
| myeloid    | C1QA              | mean           |         0     | 0.000 [0.000, 0.000]    |  0     | -0.218 [-0.341, -0.115]  | {"slide1": 0.0, "slide2": 0.0, "slide3": 0.0, "slide4": 0.0, "slide5": 0.0, "slide5b2": 0.0, "slide5b3": 0.0, "slide6": 0.0}                               |
| myeloid    | C1QB              | mean           |         0     | 0.000 [0.000, 0.000]    |  0     | -0.176 [-0.232, -0.109]  | {"slide1": 0.0, "slide2": 0.0, "slide3": 0.0, "slide4": 0.0, "slide5": 0.0, "slide5b2": 0.0, "slide5b3": 0.0, "slide6": 0.0}                               |
| myeloid    | C1QC              | mean           |         0     | 0.000 [0.000, 0.000]    |  0     | -0.353 [-0.497, -0.227]  | {"slide1": 0.0, "slide2": 0.0, "slide3": 0.0, "slide4": 0.0, "slide5": 0.0, "slide5b2": 0.0, "slide5b3": 0.0, "slide6": 0.0}                               |
| myeloid    | CD163             | mean           |         0     | 0.000 [0.000, 0.000]    |  0     | -0.211 [-0.295, -0.135]  | {"slide1": 0.0, "slide2": 0.0, "slide3": 0.0, "slide4": 0.0, "slide5": 0.0, "slide5b2": 0.0, "slide5b3": 0.0, "slide6": 0.0}                               |
| myeloid    | MRC1              | mean           |         0     | 0.000 [0.000, 0.000]    |  0     | -0.150 [-0.327, -0.053]  | {"slide1": 0.0, "slide2": 0.0, "slide3": 0.0, "slide4": 0.0, "slide5": 0.0, "slide5b2": 0.0, "slide5b3": 0.0, "slide6": 0.0}                               |
| myeloid    | SELENOP           | mean           |         0     | 0.000 [0.000, 0.000]    |  0     | -0.182 [-0.192, -0.168]  | {"slide1": 0.0, "slide2": 0.0, "slide3": 0.0, "slide4": 0.0, "slide5": 0.0, "slide5b2": 0.0, "slide5b3": 0.0, "slide6": 0.0}                               |
| myeloid    | SPP1              | mean           |         0     | 0.000 [0.000, 0.000]    |  0     | -0.549 [-0.664, -0.409]  | {"slide1": 0.0, "slide2": 0.0, "slide3": 0.0, "slide4": 0.0, "slide5": 0.0, "slide5b2": 0.0, "slide5b3": 0.0, "slide6": 0.0}                               |
| myeloid    | GPNMB             | mean           |         0     | 0.000 [0.000, 0.000]    |  0     | -0.290 [-0.420, -0.174]  | {"slide1": 0.0, "slide2": 0.0, "slide3": 0.0, "slide4": 0.0, "slide5": 0.0, "slide5b2": 0.0, "slide5b3": 0.0, "slide6": 0.0}                               |
| myeloid    | CD9               | mean           |         0     | 0.000 [0.000, 0.000]    |  0     | -0.012 [-0.095, 0.064]   | {"slide1": 0.0, "slide2": 0.0, "slide3": 0.0, "slide4": 0.0, "slide5": 0.0, "slide5b2": 0.0, "slide5b3": 0.0, "slide6": 0.0}                               |
| myeloid    | FABP5             | mean           |         0     | 0.000 [0.000, 0.000]    |  0     | -0.112 [-0.151, -0.076]  | {"slide1": 0.0, "slide2": 0.0, "slide3": 0.0, "slide4": 0.0, "slide5": 0.0, "slide5b2": 0.0, "slide5b3": 0.0, "slide6": 0.0}                               |
| myeloid    | LGALS3            | mean           |         0     | 0.000 [0.000, 0.000]    |  0     | -0.074 [-0.090, -0.050]  | {"slide1": 0.0, "slide2": 0.0, "slide3": 0.0, "slide4": 0.0, "slide5": 0.0, "slide5b2": 0.0, "slide5b3": 0.0, "slide6": 0.0}                               |
| myeloid    | genes_pooled      | mean           |         0     | 0.000 [0.000, 0.000]    |  0     | -0.251 [-0.322, -0.192]  | {"slide1": 0.0, "slide2": 0.0, "slide3": 0.0, "slide4": 0.0, "slide5": 0.0, "slide5b2": 0.0, "slide5b3": 0.0, "slide6": 0.0}                               |
| myeloid    | myeloid:C1Q_TAM   | ridge/own      |         0.299 | 0.310 [0.211, 0.463]    | -0.01  |                          | {"slide1": 0.3152, "slide2": 0.6014, "slide3": 0.2177, "slide4": 0.2571, "slide5": 0.3017, "slide5b2": 0.0893, "slide5b3": 0.2236, "slide6": 0.1959}       |
| myeloid    | myeloid:SPP1_TAM  | ridge/own      |         0.335 | 0.325 [0.203, 0.438]    |  0.01  |                          | {"slide1": 0.2451, "slide2": 0.52, "slide3": 0.2504, "slide4": 0.2232, "slide5": 0.5075, "slide5b2": 0.3922, "slide5b3": 0.285, "slide6": 0.1537}          |
| myeloid    | C1QA              | ridge/own      |         0.18  | 0.218 [0.115, 0.341]    | -0.038 |                          | {"slide1": 0.1952, "slide2": 0.4391, "slide3": 0.1518, "slide4": 0.2351, "slide5": 0.1819, "slide5b2": -0.0071, "slide5b3": 0.127, "slide6": 0.1399}       |
| myeloid    | C1QB              | ridge/own      |         0.255 | 0.176 [0.109, 0.232]    |  0.079 |                          | {"slide1": 0.1572, "slide2": 0.1804, "slide3": 0.239, "slide4": 0.2706, "slide5": 0.1557, "slide5b2": -0.019, "slide5b3": 0.0181, "slide6": 0.1695}        |
| myeloid    | C1QC              | ridge/own      |         0.237 | 0.353 [0.227, 0.497]    | -0.116 |                          | {"slide1": 0.469, "slide2": 0.599, "slide3": 0.2168, "slide4": 0.2778, "slide5": 0.2465, "slide5b2": 0.1485, "slide5b3": 0.3665, "slide6": 0.217}          |
| myeloid    | CD163             | ridge/own      |         0.186 | 0.211 [0.135, 0.295]    | -0.026 |                          | {"slide1": 0.0977, "slide2": 0.3682, "slide3": 0.1166, "slide4": 0.163, "slide5": 0.1407, "slide5b2": 0.1842, "slide5b3": 0.1763, "slide6": 0.2538}        |
| myeloid    | MRC1              | ridge/own      |         0.094 | 0.150 [0.053, 0.327]    | -0.056 |                          | {"slide1": 0.1052, "slide2": 0.4865, "slide3": 0.0576, "slide4": 0.173, "slide5": 0.0706, "slide5b2": -0.0364, "slide5b3": 0.0587, "slide6": 0.1499}       |
| myeloid    | SELENOP           | ridge/own      |         0.229 | 0.182 [0.168, 0.192]    |  0.047 |                          | {"slide1": 0.1402, "slide2": 0.1732, "slide3": 0.204, "slide4": 0.1935, "slide5": 0.2813, "slide5b2": 0.0874, "slide5b3": 0.0531, "slide6": 0.1704}        |
| myeloid    | SPP1              | ridge/own      |         0.613 | 0.549 [0.409, 0.664]    |  0.064 |                          | {"slide1": 0.4306, "slide2": 0.689, "slide3": 0.3429, "slide4": 0.6986, "slide5": 0.6207, "slide5b2": 0.6495, "slide5b3": 0.6843, "slide6": 0.3178}        |
| myeloid    | GPNMB             | ridge/own      |         0.193 | 0.290 [0.174, 0.420]    | -0.097 |                          | {"slide1": 0.3798, "slide2": 0.5132, "slide3": 0.1203, "slide4": 0.147, "slide5": 0.3474, "slide5b2": 0.3706, "slide5b3": 0.4557, "slide6": 0.1928}        |
| myeloid    | CD9               | ridge/own      |         0.106 | 0.012 [-0.064, 0.095]   |  0.094 |                          | {"slide1": 0.0313, "slide2": -0.1154, "slide3": 0.0459, "slide4": 0.1338, "slide5": 0.0875, "slide5b2": 0.0111, "slide5b3": 0.0079, "slide6": 0.1031}      |
| myeloid    | FABP5             | ridge/own      |         0.118 | 0.112 [0.076, 0.151]    |  0.005 |                          | {"slide1": 0.067, "slide2": 0.1916, "slide3": 0.0629, "slide4": 0.1005, "slide5": 0.1495, "slide5b2": 0.1494, "slide5b3": 0.13, "slide6": 0.0807}          |
| myeloid    | LGALS3            | ridge/own      |         0.136 | 0.074 [0.050, 0.090]    |  0.062 |                          | {"slide1": 0.0368, "slide2": 0.0934, "slide3": 0.0211, "slide4": 0.0787, "slide5": 0.0598, "slide5b2": 0.1422, "slide5b3": 0.0458, "slide6": 0.0881}       |
| myeloid    | genes_pooled      | ridge/own      |         0.265 | 0.251 [0.192, 0.322]    |  0.014 |                          | {"slide1": 0.2379, "slide2": 0.3851, "slide3": 0.1643, "slide4": 0.2629, "slide5": 0.2386, "slide5b2": 0.1739, "slide5b3": 0.243, "slide6": 0.1856}        |
| myeloid    | myeloid:C1Q_TAM   | ridge/genes    |         0.311 | 0.330 [0.244, 0.465]    | -0.019 | 0.020 [-0.003, 0.044]    | {"slide1": 0.2969, "slide2": 0.5984, "slide3": 0.2468, "slide4": 0.2523, "slide5": 0.3913, "slide5b2": 0.1275, "slide5b3": 0.2818, "slide6": 0.2235}       |
| myeloid    | myeloid:SPP1_TAM  | ridge/genes    |         0.349 | 0.347 [0.232, 0.445]    |  0.002 | 0.021 [-0.011, 0.052]    | {"slide1": 0.2367, "slide2": 0.4938, "slide3": 0.3085, "slide4": 0.2095, "slide5": 0.5375, "slide5b2": 0.4674, "slide5b3": 0.3323, "slide6": 0.216}        |
| myeloid    | C1QA              | ridge/genes    |         0.161 | 0.206 [0.111, 0.329]    | -0.045 | -0.012 [-0.025, -0.002]  | {"slide1": 0.1578, "slide2": 0.4344, "slide3": 0.1407, "slide4": 0.2046, "slide5": 0.1912, "slide5b2": -0.0029, "slide5b3": 0.123, "slide6": 0.1349}       |
| myeloid    | C1QB              | ridge/genes    |         0.251 | 0.180 [0.105, 0.241]    |  0.071 | 0.004 [-0.008, 0.017]    | {"slide1": 0.1506, "slide2": 0.2059, "slide3": 0.2345, "slide4": 0.2811, "slide5": 0.1829, "slide5b2": -0.0758, "slide5b3": 0.061, "slide6": 0.1772}       |
| myeloid    | C1QC              | ridge/genes    |         0.235 | 0.340 [0.221, 0.479]    | -0.105 | -0.014 [-0.031, 0.001]   | {"slide1": 0.424, "slide2": 0.5861, "slide3": 0.1969, "slide4": 0.2405, "slide5": 0.2438, "slide5b2": 0.157, "slide5b3": 0.3829, "slide6": 0.2335}         |
| myeloid    | CD163             | ridge/genes    |         0.157 | 0.181 [0.092, 0.280]    | -0.024 | -0.031 [-0.058, -0.010]  | {"slide1": 0.0859, "slide2": 0.3511, "slide3": 0.0777, "slide4": 0.0725, "slide5": 0.155, "slide5b2": 0.1572, "slide5b3": 0.137, "slide6": 0.2548}         |
| myeloid    | MRC1              | ridge/genes    |         0.056 | 0.181 [0.115, 0.316]    | -0.125 | 0.030 [-0.017, 0.066]    | {"slide1": 0.164, "slide2": 0.4844, "slide3": 0.1162, "slide4": 0.1327, "slide5": 0.163, "slide5b2": 0.0538, "slide5b3": 0.131, "slide6": 0.1295}          |
| myeloid    | SELENOP           | ridge/genes    |         0.235 | 0.187 [0.145, 0.213]    |  0.048 | 0.006 [-0.036, 0.027]    | {"slide1": 0.1414, "slide2": 0.1629, "slide3": 0.0887, "slide4": 0.2218, "slide5": 0.3335, "slide5b2": 0.0882, "slide5b3": 0.0859, "slide6": 0.1994}       |
| myeloid    | SPP1              | ridge/genes    |         0.617 | 0.574 [0.446, 0.676]    |  0.043 | 0.025 [-0.003, 0.056]    | {"slide1": 0.4956, "slide2": 0.671, "slide3": 0.3841, "slide4": 0.6882, "slide5": 0.6649, "slide5b2": 0.7277, "slide5b3": 0.7295, "slide6": 0.3588}        |
| myeloid    | GPNMB             | ridge/genes    |         0.191 | 0.304 [0.198, 0.431]    | -0.113 | 0.014 [-0.003, 0.034]    | {"slide1": 0.3472, "slide2": 0.5176, "slide3": 0.1855, "slide4": 0.1543, "slide5": 0.3588, "slide5b2": 0.3728, "slide5b3": 0.4682, "slide6": 0.2197}       |
| myeloid    | CD9               | ridge/genes    |         0.087 | 0.022 [-0.041, 0.079]   |  0.066 | 0.009 [-0.020, 0.030]    | {"slide1": 0.0468, "slide2": -0.0824, "slide3": 0.0443, "slide4": 0.0795, "slide5": 0.0651, "slide5b2": 0.0753, "slide5b3": 0.0887, "slide6": 0.0901}      |
| myeloid    | FABP5             | ridge/genes    |         0.094 | 0.095 [0.055, 0.134]    | -0.001 | -0.017 [-0.029, -0.007]  | {"slide1": 0.0201, "slide2": 0.1614, "slide3": 0.0419, "slide4": 0.0886, "slide5": 0.1485, "slide5b2": 0.1521, "slide5b3": 0.1522, "slide6": 0.0727}       |
| myeloid    | LGALS3            | ridge/genes    |         0.106 | 0.096 [0.068, 0.123]    |  0.01  | 0.022 [0.000, 0.044]     | {"slide1": 0.0628, "slide2": 0.0812, "slide3": 0.0534, "slide4": 0.0737, "slide5": 0.1072, "slide5b2": 0.148, "slide5b3": 0.1345, "slide6": 0.1465}        |
| myeloid    | genes_pooled      | ridge/genes    |         0.253 | 0.254 [0.200, 0.317]    | -0     | 0.002 [-0.009, 0.014]    | {"slide1": 0.2371, "slide2": 0.3815, "slide3": 0.1696, "slide4": 0.2424, "slide5": 0.2672, "slide5b2": 0.187, "slide5b3": 0.2692, "slide6": 0.2015}        |
| myeloid    | myeloid:C1Q_TAM   | xgb/genes      |         0.334 | 0.382 [0.304, 0.500]    | -0.048 | 0.072 [0.031, 0.107]     | {"slide1": 0.3971, "slide2": 0.6054, "slide3": 0.3251, "slide4": 0.2943, "slide5": 0.414, "slide5b2": 0.1595, "slide5b3": 0.3181, "slide6": 0.329}         |
| myeloid    | myeloid:SPP1_TAM  | xgb/genes      |         0.385 | 0.375 [0.270, 0.467]    |  0.01  | 0.050 [0.020, 0.076]     | {"slide1": 0.2502, "slide2": 0.5214, "slide3": 0.3144, "slide4": 0.2831, "slide5": 0.5416, "slide5b2": 0.4765, "slide5b3": 0.3706, "slide6": 0.2519}       |
| myeloid    | C1QA              | xgb/genes      |         0     | 0.000 [0.000, 0.000]    |  0     | -0.218 [-0.341, -0.115]  | {"slide1": 0.0, "slide2": 0.0, "slide3": 0.0, "slide4": 0.0, "slide5": 0.0, "slide5b2": 0.0, "slide5b3": 0.0, "slide6": 0.0}                               |
| myeloid    | C1QB              | xgb/genes      |         0     | 0.000 [0.000, 0.000]    |  0     | -0.176 [-0.232, -0.109]  | {"slide1": 0.0, "slide2": 0.0, "slide3": 0.0, "slide4": 0.0, "slide5": 0.0, "slide5b2": 0.0, "slide5b3": 0.0, "slide6": 0.0}                               |
| myeloid    | C1QC              | xgb/genes      |         0     | 0.000 [0.000, 0.000]    |  0     | -0.353 [-0.497, -0.227]  | {"slide1": 0.0, "slide2": 0.0, "slide3": 0.0, "slide4": 0.0, "slide5": 0.0, "slide5b2": 0.0, "slide5b3": 0.0, "slide6": 0.0}                               |
| myeloid    | CD163             | xgb/genes      |         0     | 0.000 [0.000, 0.000]    |  0     | -0.211 [-0.295, -0.135]  | {"slide1": 0.0, "slide2": 0.0, "slide3": 0.0, "slide4": 0.0, "slide5": 0.0, "slide5b2": 0.0, "slide5b3": 0.0, "slide6": 0.0}                               |
| myeloid    | MRC1              | xgb/genes      |         0     | 0.000 [0.000, 0.000]    |  0     | -0.150 [-0.327, -0.053]  | {"slide1": 0.0, "slide2": 0.0, "slide3": 0.0, "slide4": 0.0, "slide5": 0.0, "slide5b2": 0.0, "slide5b3": 0.0, "slide6": 0.0}                               |
| myeloid    | SELENOP           | xgb/genes      |         0     | 0.000 [0.000, 0.000]    |  0     | -0.182 [-0.192, -0.168]  | {"slide1": 0.0, "slide2": 0.0, "slide3": 0.0, "slide4": 0.0, "slide5": 0.0, "slide5b2": 0.0, "slide5b3": 0.0, "slide6": 0.0}                               |
| myeloid    | SPP1              | xgb/genes      |         0     | 0.000 [0.000, 0.000]    |  0     | -0.549 [-0.664, -0.409]  | {"slide1": 0.0, "slide2": 0.0, "slide3": 0.0, "slide4": 0.0, "slide5": 0.0, "slide5b2": 0.0, "slide5b3": 0.0, "slide6": 0.0}                               |
| myeloid    | GPNMB             | xgb/genes      |         0     | 0.000 [0.000, 0.000]    |  0     | -0.290 [-0.420, -0.174]  | {"slide1": 0.0, "slide2": 0.0, "slide3": 0.0, "slide4": 0.0, "slide5": 0.0, "slide5b2": 0.0, "slide5b3": 0.0, "slide6": 0.0}                               |
| myeloid    | CD9               | xgb/genes      |         0     | 0.000 [0.000, 0.000]    |  0     | -0.012 [-0.095, 0.064]   | {"slide1": 0.0, "slide2": 0.0, "slide3": 0.0, "slide4": 0.0, "slide5": 0.0, "slide5b2": 0.0, "slide5b3": 0.0, "slide6": 0.0}                               |
| myeloid    | FABP5             | xgb/genes      |         0     | 0.000 [0.000, 0.000]    |  0     | -0.112 [-0.151, -0.076]  | {"slide1": 0.0, "slide2": 0.0, "slide3": 0.0, "slide4": 0.0, "slide5": 0.0, "slide5b2": 0.0, "slide5b3": 0.0, "slide6": 0.0}                               |
| myeloid    | LGALS3            | xgb/genes      |         0     | 0.000 [0.000, 0.000]    |  0     | -0.074 [-0.090, -0.050]  | {"slide1": 0.0, "slide2": 0.0, "slide3": 0.0, "slide4": 0.0, "slide5": 0.0, "slide5b2": 0.0, "slide5b3": 0.0, "slide6": 0.0}                               |
| myeloid    | genes_pooled      | xgb/genes      |         0     | 0.000 [0.000, 0.000]    |  0     | -0.251 [-0.322, -0.192]  | {"slide1": 0.0, "slide2": 0.0, "slide3": 0.0, "slide4": 0.0, "slide5": 0.0, "slide5b2": 0.0, "slide5b3": 0.0, "slide6": 0.0}                               |
| myeloid    | myeloid:C1Q_TAM   | ridge/programs |         0.333 | 0.359 [0.271, 0.496]    | -0.026 | 0.049 [0.029, 0.068]     | {"slide1": 0.3428, "slide2": 0.6269, "slide3": 0.2477, "slide4": 0.2975, "slide5": 0.3904, "slide5b2": 0.1649, "slide5b3": 0.2948, "slide6": 0.2704}       |
| myeloid    | myeloid:SPP1_TAM  | ridge/programs |         0.383 | 0.357 [0.254, 0.442]    |  0.026 | 0.032 [-0.006, 0.062]    | {"slide1": 0.2338, "slide2": 0.4845, "slide3": 0.3138, "slide4": 0.2738, "slide5": 0.5516, "slide5b2": 0.4738, "slide5b3": 0.3467, "slide6": 0.217}        |
| myeloid    | C1QA              | ridge/programs |         0.197 | 0.236 [0.140, 0.358]    | -0.039 | 0.018 [-0.002, 0.037]    | {"slide1": 0.2439, "slide2": 0.4633, "slide3": 0.1267, "slide4": 0.227, "slide5": 0.2199, "slide5b2": 0.0415, "slide5b3": 0.1762, "slide6": 0.1572}        |
| myeloid    | C1QB              | ridge/programs |         0.261 | 0.200 [0.147, 0.252]    |  0.061 | 0.024 [0.004, 0.044]     | {"slide1": 0.2061, "slide2": 0.2018, "slide3": 0.2207, "slide4": 0.2962, "slide5": 0.2099, "slide5b2": 0.037, "slide5b3": 0.0988, "slide6": 0.1706}        |
| myeloid    | C1QC              | ridge/programs |         0.251 | 0.363 [0.238, 0.504]    | -0.112 | 0.010 [0.002, 0.019]     | {"slide1": 0.4821, "slide2": 0.606, "slide3": 0.2127, "slide4": 0.2755, "slide5": 0.2653, "slide5b2": 0.1755, "slide5b3": 0.3921, "slide6": 0.2362}        |
| myeloid    | CD163             | ridge/programs |         0.19  | 0.223 [0.142, 0.314]    | -0.033 | 0.012 [-0.007, 0.029]    | {"slide1": 0.1494, "slide2": 0.3836, "slide3": 0.075, "slide4": 0.1779, "slide5": 0.1544, "slide5b2": 0.2009, "slide5b3": 0.2267, "slide6": 0.2598}        |
| myeloid    | MRC1              | ridge/programs |         0.095 | 0.178 [0.099, 0.333]    | -0.083 | 0.028 [0.001, 0.049]     | {"slide1": 0.1407, "slide2": 0.5011, "slide3": 0.0985, "slide4": 0.1697, "slide5": 0.1362, "slide5b2": 0.0284, "slide5b3": 0.1233, "slide6": 0.1334}       |
| myeloid    | SELENOP           | ridge/programs |         0.245 | 0.196 [0.168, 0.214]    |  0.049 | 0.014 [-0.005, 0.027]    | {"slide1": 0.1482, "slide2": 0.1602, "slide3": 0.1824, "slide4": 0.2168, "slide5": 0.3314, "slide5b2": 0.0927, "slide5b3": 0.0998, "slide6": 0.1975}       |
| myeloid    | SPP1              | ridge/programs |         0.633 | 0.587 [0.461, 0.681]    |  0.046 | 0.038 [0.017, 0.060]     | {"slide1": 0.4903, "slide2": 0.6951, "slide3": 0.403, "slide4": 0.7076, "slide5": 0.6672, "slide5b2": 0.724, "slide5b3": 0.7283, "slide6": 0.3766}         |
| myeloid    | GPNMB             | ridge/programs |         0.212 | 0.318 [0.205, 0.453]    | -0.106 | 0.028 [0.015, 0.036]     | {"slide1": 0.3628, "slide2": 0.5481, "slide3": 0.1582, "slide4": 0.1833, "slide5": 0.3804, "slide5b2": 0.3975, "slide5b3": 0.4755, "slide6": 0.2231}       |
| myeloid    | CD9               | ridge/programs |         0.112 | 0.052 [-0.010, 0.117]   |  0.06  | 0.039 [0.019, 0.054]     | {"slide1": 0.0671, "slide2": -0.0527, "slide3": 0.0692, "slide4": 0.1388, "slide5": 0.1173, "slide5b2": 0.0707, "slide5b3": 0.0789, "slide6": 0.127}       |
| myeloid    | FABP5             | ridge/programs |         0.124 | 0.114 [0.069, 0.156]    |  0.01  | 0.002 [-0.016, 0.012]    | {"slide1": 0.0123, "slide2": 0.1916, "slide3": 0.0689, "slide4": 0.1059, "slide5": 0.1622, "slide5b2": 0.1688, "slide5b3": 0.1697, "slide6": 0.0907}       |
| myeloid    | LGALS3            | ridge/programs |         0.142 | 0.103 [0.069, 0.129]    |  0.039 | 0.029 [0.014, 0.046]     | {"slide1": 0.1058, "slide2": 0.1448, "slide3": 0.0201, "slide4": 0.0961, "slide5": 0.08, "slide5b2": 0.1714, "slide5b3": 0.1477, "slide6": 0.1017}         |
| myeloid    | genes_pooled      | ridge/programs |         0.277 | 0.274 [0.212, 0.341]    |  0.003 | 0.022 [0.013, 0.033]     | {"slide1": 0.2717, "slide2": 0.4051, "slide3": 0.1704, "slide4": 0.2751, "slide5": 0.2759, "slide5b2": 0.2161, "slide5b3": 0.2941, "slide6": 0.206}        |
| myeloid    | myeloid:C1Q_TAM   | xgb/programs   |         0.331 | 0.391 [0.308, 0.520]    | -0.06  | 0.081 [0.051, 0.111]     | {"slide1": 0.3968, "slide2": 0.6375, "slide3": 0.3337, "slide4": 0.2975, "slide5": 0.407, "slide5b2": 0.1774, "slide5b3": 0.3246, "slide6": 0.3276}        |
| myeloid    | myeloid:SPP1_TAM  | xgb/programs   |         0.387 | 0.360 [0.257, 0.444]    |  0.027 | 0.035 [-0.002, 0.067]    | {"slide1": 0.2201, "slide2": 0.4918, "slide3": 0.3101, "slide4": 0.2803, "slide5": 0.5414, "slide5b2": 0.4702, "slide5b3": 0.3307, "slide6": 0.2377}       |
| myeloid    | C1QA              | xgb/programs   |         0     | 0.000 [0.000, 0.000]    |  0     | -0.218 [-0.341, -0.115]  | {"slide1": 0.0, "slide2": 0.0, "slide3": 0.0, "slide4": 0.0, "slide5": 0.0, "slide5b2": 0.0, "slide5b3": 0.0, "slide6": 0.0}                               |
| myeloid    | C1QB              | xgb/programs   |         0     | 0.000 [0.000, 0.000]    |  0     | -0.176 [-0.232, -0.109]  | {"slide1": 0.0, "slide2": 0.0, "slide3": 0.0, "slide4": 0.0, "slide5": 0.0, "slide5b2": 0.0, "slide5b3": 0.0, "slide6": 0.0}                               |
| myeloid    | C1QC              | xgb/programs   |         0     | 0.000 [0.000, 0.000]    |  0     | -0.353 [-0.497, -0.227]  | {"slide1": 0.0, "slide2": 0.0, "slide3": 0.0, "slide4": 0.0, "slide5": 0.0, "slide5b2": 0.0, "slide5b3": 0.0, "slide6": 0.0}                               |
| myeloid    | CD163             | xgb/programs   |         0     | 0.000 [0.000, 0.000]    |  0     | -0.211 [-0.295, -0.135]  | {"slide1": 0.0, "slide2": 0.0, "slide3": 0.0, "slide4": 0.0, "slide5": 0.0, "slide5b2": 0.0, "slide5b3": 0.0, "slide6": 0.0}                               |
| myeloid    | MRC1              | xgb/programs   |         0     | 0.000 [0.000, 0.000]    |  0     | -0.150 [-0.327, -0.053]  | {"slide1": 0.0, "slide2": 0.0, "slide3": 0.0, "slide4": 0.0, "slide5": 0.0, "slide5b2": 0.0, "slide5b3": 0.0, "slide6": 0.0}                               |
| myeloid    | SELENOP           | xgb/programs   |         0     | 0.000 [0.000, 0.000]    |  0     | -0.182 [-0.192, -0.168]  | {"slide1": 0.0, "slide2": 0.0, "slide3": 0.0, "slide4": 0.0, "slide5": 0.0, "slide5b2": 0.0, "slide5b3": 0.0, "slide6": 0.0}                               |
| myeloid    | SPP1              | xgb/programs   |         0     | 0.000 [0.000, 0.000]    |  0     | -0.549 [-0.664, -0.409]  | {"slide1": 0.0, "slide2": 0.0, "slide3": 0.0, "slide4": 0.0, "slide5": 0.0, "slide5b2": 0.0, "slide5b3": 0.0, "slide6": 0.0}                               |
| myeloid    | GPNMB             | xgb/programs   |         0     | 0.000 [0.000, 0.000]    |  0     | -0.290 [-0.420, -0.174]  | {"slide1": 0.0, "slide2": 0.0, "slide3": 0.0, "slide4": 0.0, "slide5": 0.0, "slide5b2": 0.0, "slide5b3": 0.0, "slide6": 0.0}                               |
| myeloid    | CD9               | xgb/programs   |         0     | 0.000 [0.000, 0.000]    |  0     | -0.012 [-0.095, 0.064]   | {"slide1": 0.0, "slide2": 0.0, "slide3": 0.0, "slide4": 0.0, "slide5": 0.0, "slide5b2": 0.0, "slide5b3": 0.0, "slide6": 0.0}                               |
| myeloid    | FABP5             | xgb/programs   |         0     | 0.000 [0.000, 0.000]    |  0     | -0.112 [-0.151, -0.076]  | {"slide1": 0.0, "slide2": 0.0, "slide3": 0.0, "slide4": 0.0, "slide5": 0.0, "slide5b2": 0.0, "slide5b3": 0.0, "slide6": 0.0}                               |
| myeloid    | LGALS3            | xgb/programs   |         0     | 0.000 [0.000, 0.000]    |  0     | -0.074 [-0.090, -0.050]  | {"slide1": 0.0, "slide2": 0.0, "slide3": 0.0, "slide4": 0.0, "slide5": 0.0, "slide5b2": 0.0, "slide5b3": 0.0, "slide6": 0.0}                               |
| myeloid    | genes_pooled      | xgb/programs   |         0     | 0.000 [0.000, 0.000]    |  0     | -0.251 [-0.322, -0.192]  | {"slide1": 0.0, "slide2": 0.0, "slide3": 0.0, "slide4": 0.0, "slide5": 0.0, "slide5b2": 0.0, "slide5b3": 0.0, "slide6": 0.0}                               |
| myeloid    | myeloid:C1Q_TAM   | graph/sage     |         0.336 | 0.324 [0.239, 0.452]    |  0.013 | 0.014 [-0.014, 0.046]    | {"slide1": 0.3213, "slide2": 0.5656, "slide3": 0.2238, "slide4": 0.2725, "slide5": 0.3496, "slide5b2": 0.0626, "slide5b3": 0.311, "slide6": 0.2763}        |
| myeloid    | myeloid:SPP1_TAM  | graph/sage     |         0.377 | 0.375 [0.272, 0.466]    |  0.002 | 0.050 [0.016, 0.086]     | {"slide1": 0.2921, "slide2": 0.5187, "slide3": 0.32, "slide4": 0.239, "slide5": 0.5314, "slide5b2": 0.4707, "slide5b3": 0.3736, "slide6": 0.2742}          |
| myeloid    | C1QA              | graph/sage     |         0.192 | 0.202 [0.103, 0.339]    | -0.01  | -0.016 [-0.043, 0.005]   | {"slide1": 0.2127, "slide2": 0.4503, "slide3": 0.1196, "slide4": 0.1659, "slide5": 0.1889, "slide5b2": -0.028, "slide5b3": 0.1501, "slide6": 0.1239}       |
| myeloid    | C1QB              | graph/sage     |         0.279 | 0.156 [0.091, 0.214]    |  0.123 | -0.020 [-0.062, 0.024]   | {"slide1": 0.1576, "slide2": 0.242, "slide3": 0.1324, "slide4": 0.2088, "slide5": 0.165, "slide5b2": -0.0869, "slide5b3": 0.0925, "slide6": 0.1626}        |
| myeloid    | C1QC              | graph/sage     |         0.265 | 0.321 [0.193, 0.464]    | -0.057 | -0.032 [-0.051, -0.014]  | {"slide1": 0.4506, "slide2": 0.5687, "slide3": 0.166, "slide4": 0.2052, "slide5": 0.2361, "slide5b2": 0.1109, "slide5b3": 0.3782, "slide6": 0.2224}        |
| myeloid    | CD163             | graph/sage     |         0.195 | 0.219 [0.139, 0.310]    | -0.024 | 0.008 [-0.007, 0.021]    | {"slide1": 0.1193, "slide2": 0.3778, "slide3": 0.109, "slide4": 0.1417, "slide5": 0.1626, "slide5b2": 0.1958, "slide5b3": 0.2331, "slide6": 0.2824}        |
| myeloid    | MRC1              | graph/sage     |         0.092 | 0.182 [0.103, 0.329]    | -0.09  | 0.032 [0.004, 0.055]     | {"slide1": 0.1751, "slide2": 0.4853, "slide3": 0.0888, "slide4": 0.1617, "slide5": 0.1238, "slide5b2": 0.0349, "slide5b3": 0.1317, "slide6": 0.1785}       |
| myeloid    | SELENOP           | graph/sage     |         0.271 | 0.227 [0.200, 0.246]    |  0.044 | 0.045 [0.023, 0.059]     | {"slide1": 0.1712, "slide2": 0.2017, "slide3": 0.2011, "slide4": 0.2432, "slide5": 0.3563, "slide5b2": 0.134, "slide5b3": 0.1799, "slide6": 0.2251}        |
| myeloid    | SPP1              | graph/sage     |         0.654 | 0.627 [0.521, 0.714]    |  0.027 | 0.078 [0.048, 0.110]     | {"slide1": 0.522, "slide2": 0.717, "slide3": 0.4623, "slide4": 0.7527, "slide5": 0.6991, "slide5b2": 0.7405, "slide5b3": 0.7451, "slide6": 0.4489}         |
| myeloid    | GPNMB             | graph/sage     |         0.228 | 0.295 [0.168, 0.441]    | -0.067 | 0.005 [-0.014, 0.024]    | {"slide1": 0.3985, "slide2": 0.5419, "slide3": 0.1492, "slide4": 0.121, "slide5": 0.3498, "slide5b2": 0.3745, "slide5b3": 0.4729, "slide6": 0.178}         |
| myeloid    | CD9               | graph/sage     |         0.124 | 0.038 [-0.036, 0.122]   |  0.086 | 0.026 [0.001, 0.052]     | {"slide1": 0.1195, "slide2": -0.0872, "slide3": 0.0403, "slide4": 0.1571, "slide5": 0.1124, "slide5b2": 0.084, "slide5b3": 0.1251, "slide6": 0.0819}       |
| myeloid    | FABP5             | graph/sage     |         0.141 | 0.163 [0.114, 0.217]    | -0.021 | 0.050 [0.033, 0.070]     | {"slide1": 0.1627, "slide2": 0.2636, "slide3": 0.0878, "slide4": 0.1337, "slide5": 0.1896, "slide5b2": 0.2248, "slide5b3": 0.2352, "slide6": 0.1057}       |
| myeloid    | LGALS3            | graph/sage     |         0.152 | 0.108 [0.053, 0.171]    |  0.044 | 0.034 [-0.017, 0.091]    | {"slide1": 0.138, "slide2": 0.2001, "slide3": 0.0395, "slide4": 0.0533, "slide5": 0.1034, "slide5b2": 0.2654, "slide5b3": 0.2255, "slide6": 0.04}          |
| myeloid    | genes_pooled      | graph/sage     |         0.289 | 0.266 [0.203, 0.339]    |  0.024 | 0.015 [-0.002, 0.029]    | {"slide1": 0.2766, "slide2": 0.4098, "slide3": 0.1636, "slide4": 0.2433, "slide5": 0.271, "slide5b2": 0.192, "slide5b3": 0.3047, "slide6": 0.2038}         |
| tumour     | tumour:classical  | mean           |         0     | 0.000 [0.000, 0.000]    |  0     | 0.074 [-0.266, 0.346]    | {"slide1": 0.0, "slide2": 0.0, "slide3": 0.0, "slide4": 0.0, "slide5": 0.0, "slide5b2": 0.0, "slide5b3": 0.0, "slide6": 0.0}                               |
| tumour     | tumour:basal_like | mean           |         0     | 0.000 [0.000, 0.000]    |  0     | -0.375 [-0.433, -0.290]  | {"slide1": 0.0, "slide2": 0.0, "slide3": 0.0, "slide4": 0.0, "slide5": 0.0, "slide5b2": 0.0, "slide5b3": 0.0, "slide6": 0.0}                               |
| tumour     | CEACAM6           | mean           |         0     | 0.000 [0.000, 0.000]    |  0     | -0.285 [-0.411, -0.192]  | {"slide1": 0.0, "slide2": 0.0, "slide3": 0.0, "slide4": 0.0, "slide5": 0.0, "slide5b2": 0.0, "slide5b3": 0.0, "slide6": 0.0}                               |
| tumour     | AGR2              | mean           |         0     | 0.000 [0.000, 0.000]    |  0     | -0.275 [-0.352, -0.224]  | {"slide1": 0.0, "slide2": 0.0, "slide3": 0.0, "slide4": 0.0, "slide5": 0.0, "slide5b2": 0.0, "slide5b3": 0.0, "slide6": 0.0}                               |
| tumour     | S100P             | mean           |         0     | 0.000 [0.000, 0.000]    |  0     | -0.323 [-0.412, -0.250]  | {"slide1": 0.0, "slide2": 0.0, "slide3": 0.0, "slide4": 0.0, "slide5": 0.0, "slide5b2": 0.0, "slide5b3": 0.0, "slide6": 0.0}                               |
| tumour     | LCN2              | mean           |         0     | 0.000 [0.000, 0.000]    |  0     | -0.414 [-0.588, -0.291]  | {"slide1": 0.0, "slide2": 0.0, "slide3": 0.0, "slide4": 0.0, "slide5": 0.0, "slide5b2": 0.0, "slide5b3": 0.0, "slide6": 0.0}                               |
| tumour     | PIGR              | mean           |         0     | 0.000 [0.000, 0.000]    |  0     | -0.270 [-0.345, -0.160]  | {"slide1": 0.0, "slide2": 0.0, "slide3": 0.0, "slide4": 0.0, "slide5": 0.0, "slide5b2": 0.0, "slide5b3": 0.0, "slide6": 0.0}                               |
| tumour     | OLFM4             | mean           |         0     | 0.000 [0.000, 0.000]    |  0     | -0.244 [-0.291, -0.132]  | {"slide1": 0.0, "slide2": 0.0, "slide3": 0.0, "slide4": 0.0, "slide5": 0.0, "slide5b2": 0.0, "slide5b3": 0.0, "slide6": 0.0}                               |
| tumour     | SPINK1            | mean           |         0     | 0.000 [0.000, 0.000]    |  0     | -0.458 [-0.514, -0.331]  | {"slide1": 0.0, "slide2": 0.0, "slide3": 0.0, "slide4": 0.0, "slide5": 0.0, "slide5b2": 0.0, "slide5b3": 0.0, "slide6": 0.0}                               |
| tumour     | CEACAM1           | mean           |         0     | 0.000 [0.000, 0.000]    |  0     | 0.005 [-0.043, 0.052]    | {"slide1": 0.0, "slide2": 0.0, "slide3": 0.0, "slide4": 0.0, "slide5": 0.0, "slide5b2": 0.0, "slide5b3": 0.0, "slide6": 0.0}                               |
| tumour     | KRT6A             | mean           |         0     | 0.000 [0.000, 0.000]    |  0     | -0.314 [-0.379, -0.233]  | {"slide1": 0.0, "slide2": 0.0, "slide3": 0.0, "slide4": 0.0, "slide5": 0.0, "slide5b2": 0.0, "slide5b3": 0.0, "slide6": 0.0}                               |
| tumour     | KRT14             | mean           |         0     | 0.000 [0.000, 0.000]    |  0     | -0.275 [-0.342, -0.183]  | {"slide1": 0.0, "slide2": 0.0, "slide3": 0.0, "slide4": 0.0, "slide5": 0.0, "slide5b2": 0.0, "slide5b3": 0.0, "slide6": 0.0}                               |
| tumour     | KRT17             | mean           |         0     | 0.000 [0.000, 0.000]    |  0     | -0.471 [-0.507, -0.411]  | {"slide1": 0.0, "slide2": 0.0, "slide3": 0.0, "slide4": 0.0, "slide5": 0.0, "slide5b2": 0.0, "slide5b3": 0.0, "slide6": 0.0}                               |
| tumour     | S100A2            | mean           |         0     | 0.000 [0.000, 0.000]    |  0     | -0.464 [-0.592, -0.280]  | {"slide1": 0.0, "slide2": 0.0, "slide3": 0.0, "slide4": 0.0, "slide5": 0.0, "slide5b2": 0.0, "slide5b3": 0.0, "slide6": 0.0}                               |
| tumour     | SERPINB5          | mean           |         0     | 0.000 [0.000, 0.000]    |  0     | -0.183 [-0.236, -0.125]  | {"slide1": 0.0, "slide2": 0.0, "slide3": 0.0, "slide4": 0.0, "slide5": 0.0, "slide5b2": 0.0, "slide5b3": 0.0, "slide6": 0.0}                               |
| tumour     | AREG              | mean           |         0     | 0.000 [0.000, 0.000]    |  0     | -0.137 [-0.184, -0.101]  | {"slide1": 0.0, "slide2": 0.0, "slide3": 0.0, "slide4": 0.0, "slide5": 0.0, "slide5b2": 0.0, "slide5b3": 0.0, "slide6": 0.0}                               |
| tumour     | SLC2A1            | mean           |         0     | 0.000 [0.000, 0.000]    |  0     | -0.159 [-0.223, -0.120]  | {"slide1": 0.0, "slide2": 0.0, "slide3": 0.0, "slide4": 0.0, "slide5": 0.0, "slide5b2": 0.0, "slide5b3": 0.0, "slide6": 0.0}                               |
| tumour     | KRT15             | mean           |         0     | 0.000 [0.000, 0.000]    |  0     | 0.022 [-0.006, 0.057]    | {"slide1": 0.0, "slide2": 0.0, "slide3": 0.0, "slide4": 0.0, "slide5": 0.0, "slide5b2": 0.0, "slide5b3": 0.0, "slide6": 0.0}                               |
| tumour     | ITGA3             | mean           |         0     | 0.000 [0.000, 0.000]    |  0     | -0.074 [-0.167, 0.001]   | {"slide1": 0.0, "slide2": 0.0, "slide3": 0.0, "slide4": 0.0, "slide5": 0.0, "slide5b2": 0.0, "slide5b3": 0.0, "slide6": 0.0}                               |
| tumour     | MET               | mean           |         0     | 0.000 [0.000, 0.000]    |  0     | -0.117 [-0.154, -0.089]  | {"slide1": 0.0, "slide2": 0.0, "slide3": 0.0, "slide4": 0.0, "slide5": 0.0, "slide5b2": 0.0, "slide5b3": 0.0, "slide6": 0.0}                               |
| tumour     | genes_pooled      | mean           |         0     | 0.000 [0.000, 0.000]    |  0     | -0.286 [-0.327, -0.259]  | {"slide1": 0.0, "slide2": 0.0, "slide3": 0.0, "slide4": 0.0, "slide5": 0.0, "slide5b2": 0.0, "slide5b3": 0.0, "slide6": 0.0}                               |
| tumour     | tumour:classical  | ridge/own      |         0.618 | -0.074 [-0.346, 0.266]  |  0.692 |                          | {"slide1": 0.4658, "slide2": 0.2924, "slide3": -0.1959, "slide4": -0.3676, "slide5": -0.5256, "slide5b2": -0.5526, "slide5b3": -0.2652, "slide6": -0.2626} |
| tumour     | tumour:basal_like | ridge/own      |         0.397 | 0.375 [0.290, 0.433]    |  0.022 |                          | {"slide1": 0.4345, "slide2": 0.4959, "slide3": 0.3838, "slide4": 0.2253, "slide5": 0.4623, "slide5b2": 0.3251, "slide5b3": 0.3755, "slide6": 0.1833}       |
| tumour     | CEACAM6           | ridge/own      |         0.627 | 0.285 [0.192, 0.411]    |  0.342 |                          | {"slide1": 0.5384, "slide2": 0.1781, "slide3": 0.1063, "slide4": 0.3587, "slide5": 0.1601, "slide5b2": -0.014, "slide5b3": 0.3868, "slide6": 0.2954}       |
| tumour     | AGR2              | ridge/own      |         0.451 | 0.275 [0.224, 0.352]    |  0.176 |                          | {"slide1": 0.2815, "slide2": 0.332, "slide3": 0.3768, "slide4": 0.3562, "slide5": 0.5032, "slide5b2": 0.1271, "slide5b3": -0.0197, "slide6": 0.124}        |
| tumour     | S100P             | ridge/own      |         0.466 | 0.323 [0.250, 0.412]    |  0.143 |                          | {"slide1": 0.5415, "slide2": 0.2846, "slide3": 0.2382, "slide4": 0.2278, "slide5": 0.5017, "slide5b2": 0.3428, "slide5b3": 0.154, "slide6": 0.2749}        |
| tumour     | LCN2              | ridge/own      |         0.54  | 0.414 [0.291, 0.588]    |  0.126 |                          | {"slide1": 0.5202, "slide2": 0.4794, "slide3": 0.2864, "slide4": 0.7375, "slide5": 0.2035, "slide5b2": 0.2123, "slide5b3": 0.2799, "slide6": 0.4282}       |
| tumour     | PIGR              | ridge/own      |         0.535 | 0.270 [0.160, 0.345]    |  0.265 |                          | {"slide1": 0.2463, "slide2": 0.0902, "slide3": 0.4266, "slide4": 0.1372, "slide5": 0.4815, "slide5b2": 0.2936, "slide5b3": 0.0837, "slide6": 0.2629}       |
| tumour     | OLFM4             | ridge/own      |         0.604 | 0.244 [0.132, 0.291]    |  0.36  |                          | {"slide1": 0.1688, "slide2": 0.1406, "slide3": 0.0611, "slide4": 0.3084, "slide5": 0.4714, "slide5b2": 0.2915, "slide5b3": 0.1234, "slide6": 0.287}        |
| tumour     | SPINK1            | ridge/own      |         0.579 | 0.458 [0.331, 0.514]    |  0.121 |                          | {"slide1": 0.5589, "slide2": 0.153, "slide3": 0.5318, "slide4": 0.4459, "slide5": 0.6892, "slide5b2": 0.4064, "slide5b3": 0.3815, "slide6": 0.4855}        |
| tumour     | CEACAM1           | ridge/own      |         0.132 | -0.005 [-0.052, 0.043]  |  0.138 |                          | {"slide1": 0.0499, "slide2": -0.0006, "slide3": -0.0982, "slide4": -0.0611, "slide5": 0.0862, "slide5b2": -0.0091, "slide5b3": -0.1272, "slide6": 0.0868}  |
| tumour     | KRT6A             | ridge/own      |         0.574 | 0.314 [0.233, 0.379]    |  0.259 |                          | {"slide1": 0.3899, "slide2": 0.1483, "slide3": 0.425, "slide4": 0.2815, "slide5": 0.3937, "slide5b2": 0.2319, "slide5b3": 0.2031, "slide6": 0.2838}        |
| tumour     | KRT14             | ridge/own      |         0.307 | 0.275 [0.183, 0.342]    |  0.032 |                          | {"slide1": 0.1809, "slide2": 0.2067, "slide3": 0.1377, "slide4": 0.4004, "slide5": 0.4365, "slide5b2": 0.3205, "slide5b3": 0.3287, "slide6": 0.2539}       |
| tumour     | KRT17             | ridge/own      |         0.428 | 0.471 [0.411, 0.507]    | -0.043 |                          | {"slide1": 0.5442, "slide2": 0.5138, "slide3": 0.2532, "slide4": 0.4454, "slide5": 0.4415, "slide5b2": 0.3098, "slide5b3": 0.6051, "slide6": 0.4446}       |
| tumour     | S100A2            | ridge/own      |         0.453 | 0.464 [0.280, 0.592]    | -0.011 |                          | {"slide1": 0.6555, "slide2": 0.3394, "slide3": 0.2367, "slide4": 0.4076, "slide5": 0.3436, "slide5b2": 0.3385, "slide5b3": 0.4064, "slide6": 0.2432}       |
| tumour     | SERPINB5          | ridge/own      |         0.147 | 0.183 [0.125, 0.236]    | -0.036 |                          | {"slide1": 0.1647, "slide2": 0.2514, "slide3": 0.0534, "slide4": 0.1351, "slide5": 0.113, "slide5b2": 0.2736, "slide5b3": 0.4468, "slide6": 0.1623}        |
| tumour     | AREG              | ridge/own      |         0.157 | 0.137 [0.101, 0.184]    |  0.021 |                          | {"slide1": 0.1672, "slide2": 0.256, "slide3": 0.1824, "slide4": 0.1155, "slide5": 0.0957, "slide5b2": 0.0807, "slide5b3": 0.0882, "slide6": 0.0906}        |
| tumour     | SLC2A1            | ridge/own      |         0.217 | 0.159 [0.120, 0.223]    |  0.059 |                          | {"slide1": 0.0916, "slide2": 0.3208, "slide3": 0.1968, "slide4": 0.1015, "slide5": 0.1737, "slide5b2": 0.0651, "slide5b3": 0.2005, "slide6": 0.1612}       |
| tumour     | KRT15             | ridge/own      |         0.089 | -0.022 [-0.057, 0.006]  |  0.112 |                          | {"slide1": -0.0725, "slide2": -0.0247, "slide3": -0.0189, "slide4": -0.0607, "slide5": -0.0565, "slide5b2": 0.0665, "slide5b3": 0.0528, "slide6": -0.0464} |
| tumour     | ITGA3             | ridge/own      |         0.238 | 0.074 [-0.001, 0.167]   |  0.164 |                          | {"slide1": 0.1895, "slide2": 0.2637, "slide3": -0.1041, "slide4": 0.0313, "slide5": 0.0184, "slide5b2": -0.0186, "slide5b3": 0.1241, "slide6": -0.0019}    |
| tumour     | MET               | ridge/own      |         0.086 | 0.117 [0.089, 0.154]    | -0.03  |                          | {"slide1": 0.1229, "slide2": 0.1679, "slide3": 0.0458, "slide4": 0.1338, "slide5": 0.0486, "slide5b2": 0.0699, "slide5b3": 0.2043, "slide6": 0.1632}       |
| tumour     | genes_pooled      | ridge/own      |         0.425 | 0.286 [0.259, 0.327]    |  0.139 |                          | {"slide1": 0.3834, "slide2": 0.259, "slide3": 0.2377, "slide4": 0.3025, "slide5": 0.3637, "slide5b2": 0.2107, "slide5b3": 0.2638, "slide6": 0.2585}        |
| tumour     | tumour:classical  | ridge/genes    |         0.668 | 0.468 [0.374, 0.571]    |  0.201 | 0.542 [0.311, 0.735]     | {"slide1": 0.6425, "slide2": 0.5597, "slide3": 0.5227, "slide4": 0.2679, "slide5": 0.338, "slide5b2": 0.3427, "slide5b3": 0.4071, "slide6": 0.3816}        |
| tumour     | tumour:basal_like | ridge/genes    |         0.39  | 0.261 [0.121, 0.373]    |  0.128 | -0.114 [-0.180, -0.033]  | {"slide1": 0.4383, "slide2": 0.3679, "slide3": 0.3999, "slide4": 0.0645, "slide5": 0.3026, "slide5b2": 0.2365, "slide5b3": 0.1482, "slide6": -0.1126}      |
| tumour     | CEACAM6           | ridge/genes    |         0.646 | 0.492 [0.386, 0.602]    |  0.155 | 0.207 [0.133, 0.276]     | {"slide1": 0.6698, "slide2": 0.2536, "slide3": 0.4797, "slide4": 0.5048, "slide5": 0.4458, "slide5b2": 0.2708, "slide5b3": 0.5671, "slide6": 0.6157}       |
| tumour     | AGR2              | ridge/genes    |         0.441 | 0.449 [0.394, 0.515]    | -0.008 | 0.174 [0.123, 0.197]     | {"slide1": 0.4332, "slide2": 0.3939, "slide3": 0.5896, "slide4": 0.5006, "slide5": 0.6359, "slide5b2": 0.3375, "slide5b3": 0.2424, "slide6": 0.324}        |
| tumour     | S100P             | ridge/genes    |         0.459 | 0.391 [0.289, 0.445]    |  0.068 | 0.068 [-0.007, 0.113]    | {"slide1": 0.4926, "slide2": 0.2702, "slide3": 0.3697, "slide4": 0.2322, "slide5": 0.6148, "slide5b2": 0.4408, "slide5b3": 0.2974, "slide6": 0.4127}       |
| tumour     | LCN2              | ridge/genes    |         0.547 | 0.441 [0.340, 0.597]    |  0.106 | 0.027 [-0.013, 0.058]    | {"slide1": 0.4949, "slide2": 0.4548, "slide3": 0.3887, "slide4": 0.7261, "slide5": 0.2487, "slide5b2": 0.2547, "slide5b3": 0.3574, "slide6": 0.4772}       |
| tumour     | PIGR              | ridge/genes    |         0.532 | 0.400 [0.263, 0.471]    |  0.133 | 0.130 [0.095, 0.146]     | {"slide1": 0.3256, "slide2": 0.1575, "slide3": 0.5545, "slide4": 0.2589, "slide5": 0.6509, "slide5b2": 0.4772, "slide5b3": 0.1956, "slide6": 0.4026}       |
| tumour     | OLFM4             | ridge/genes    |         0.642 | 0.393 [0.270, 0.446]    |  0.249 | 0.150 [0.110, 0.204]     | {"slide1": 0.2743, "slide2": 0.225, "slide3": 0.306, "slide4": 0.3903, "slide5": 0.5929, "slide5b2": 0.4492, "slide5b3": 0.2746, "slide6": 0.5255}         |
| tumour     | SPINK1            | ridge/genes    |         0.595 | 0.527 [0.370, 0.602]    |  0.068 | 0.069 [0.033, 0.096]     | {"slide1": 0.5991, "slide2": 0.1289, "slide3": 0.6369, "slide4": 0.5299, "slide5": 0.7226, "slide5b2": 0.4719, "slide5b3": 0.4828, "slide6": 0.6007}       |
| tumour     | CEACAM1           | ridge/genes    |         0.091 | 0.042 [-0.031, 0.097]   |  0.049 | 0.048 [-0.012, 0.098]    | {"slide1": 0.0849, "slide2": -0.0507, "slide3": 0.079, "slide4": -0.0932, "slide5": 0.1532, "slide5b2": 0.0614, "slide5b3": -0.0694, "slide6": 0.1463}     |
| tumour     | KRT6A             | ridge/genes    |         0.565 | 0.296 [0.168, 0.414]    |  0.269 | -0.019 [-0.081, 0.035]   | {"slide1": 0.4292, "slide2": 0.1357, "slide3": 0.4712, "slide4": 0.2632, "slide5": 0.32, "slide5b2": 0.0692, "slide5b3": 0.068, "slide6": 0.3113}          |
| tumour     | KRT14             | ridge/genes    |         0.255 | 0.060 [-0.048, 0.137]   |  0.195 | -0.215 [-0.308, -0.167]  | {"slide1": 0.0428, "slide2": -0.0699, "slide3": -0.0005, "slide4": 0.0592, "slide5": 0.3941, "slide5b2": 0.1504, "slide5b3": 0.0392, "slide6": -0.1021}    |
| tumour     | KRT17             | ridge/genes    |         0.413 | 0.451 [0.392, 0.479]    | -0.038 | -0.020 [-0.091, 0.026]   | {"slide1": 0.3674, "slide2": 0.4684, "slide3": 0.3543, "slide4": 0.4537, "slide5": 0.4666, "slide5b2": 0.3267, "slide5b3": 0.5999, "slide6": 0.4441}       |
| tumour     | S100A2            | ridge/genes    |         0.42  | 0.425 [0.223, 0.573]    | -0.006 | -0.039 [-0.107, -0.003]  | {"slide1": 0.6428, "slide2": 0.2209, "slide3": 0.2946, "slide4": 0.2405, "slide5": 0.3011, "slide5b2": 0.2809, "slide5b3": 0.3299, "slide6": 0.1813}       |
| tumour     | SERPINB5          | ridge/genes    |         0.086 | 0.088 [0.010, 0.138]    | -0.002 | -0.095 [-0.129, -0.061]  | {"slide1": 0.0999, "slide2": 0.1181, "slide3": -0.1074, "slide4": 0.0298, "slide5": 0.0872, "slide5b2": 0.1153, "slide5b3": 0.2769, "slide6": 0.1377}      |
| tumour     | AREG              | ridge/genes    |         0.085 | 0.083 [0.012, 0.142]    |  0.002 | -0.054 [-0.122, -0.021]  | {"slide1": 0.1274, "slide2": 0.0427, "slide3": 0.1908, "slide4": 0.002, "slide5": 0.0788, "slide5b2": 0.0663, "slide5b3": 0.0037, "slide6": -0.0218}       |
| tumour     | SLC2A1            | ridge/genes    |         0.191 | -0.058 [-0.109, -0.012] |  0.248 | -0.216 [-0.286, -0.156]  | {"slide1": 0.0123, "slide2": -0.0251, "slide3": -0.0814, "slide4": -0.1783, "slide5": -0.0999, "slide5b2": -0.0909, "slide5b3": -0.0019, "slide6": 0.0116} |
| tumour     | KRT15             | ridge/genes    |         0.031 | -0.051 [-0.102, -0.011] |  0.081 | -0.029 [-0.069, 0.024]   | {"slide1": -0.0615, "slide2": -0.0774, "slide3": 0.056, "slide4": -0.0281, "slide5": -0.025, "slide5b2": -0.0235, "slide5b3": -0.0538, "slide6": -0.1766}  |
| tumour     | ITGA3             | ridge/genes    |         0.201 | 0.055 [-0.067, 0.173]   |  0.146 | -0.019 [-0.124, 0.053]   | {"slide1": 0.3197, "slide2": 0.1532, "slide3": 0.0252, "slide4": -0.0994, "slide5": -0.0139, "slide5b2": -0.0213, "slide5b3": 0.173, "slide6": -0.203}     |
| tumour     | MET               | ridge/genes    |         0.036 | 0.080 [0.056, 0.108]    | -0.044 | -0.037 [-0.074, -0.016]  | {"slide1": 0.0466, "slide2": 0.0589, "slide3": 0.0588, "slide4": 0.062, "slide5": 0.0384, "slide5b2": 0.0527, "slide5b3": 0.1701, "slide6": 0.1433}        |
| tumour     | genes_pooled      | ridge/genes    |         0.412 | 0.325 [0.276, 0.352]    |  0.088 | 0.038 [-0.005, 0.066]    | {"slide1": 0.3893, "slide2": 0.2218, "slide3": 0.3496, "slide4": 0.3019, "slide5": 0.4265, "slide5b2": 0.2718, "slide5b3": 0.3226, "slide6": 0.3072}       |
| tumour     | tumour:classical  | xgb/genes      |         0.694 | 0.587 [0.493, 0.707]    |  0.107 | 0.661 [0.449, 0.842]     | {"slide1": 0.7993, "slide2": 0.6616, "slide3": 0.5665, "slide4": 0.4979, "slide5": 0.4232, "slide5b2": 0.4651, "slide5b3": 0.4764, "slide6": 0.546}        |
| tumour     | tumour:basal_like | xgb/genes      |         0.444 | 0.500 [0.419, 0.546]    | -0.056 | 0.125 [0.100, 0.143]     | {"slide1": 0.5143, "slide2": 0.592, "slide3": 0.5353, "slide4": 0.3638, "slide5": 0.5519, "slide5b2": 0.4613, "slide5b3": 0.5764, "slide6": 0.3241}        |
| tumour     | CEACAM6           | xgb/genes      |         0     | 0.000 [0.000, 0.000]    |  0     | -0.285 [-0.411, -0.192]  | {"slide1": 0.0, "slide2": 0.0, "slide3": 0.0, "slide4": 0.0, "slide5": 0.0, "slide5b2": 0.0, "slide5b3": 0.0, "slide6": 0.0}                               |
| tumour     | AGR2              | xgb/genes      |         0     | 0.000 [0.000, 0.000]    |  0     | -0.275 [-0.352, -0.224]  | {"slide1": 0.0, "slide2": 0.0, "slide3": 0.0, "slide4": 0.0, "slide5": 0.0, "slide5b2": 0.0, "slide5b3": 0.0, "slide6": 0.0}                               |
| tumour     | S100P             | xgb/genes      |         0     | 0.000 [0.000, 0.000]    |  0     | -0.323 [-0.412, -0.250]  | {"slide1": 0.0, "slide2": 0.0, "slide3": 0.0, "slide4": 0.0, "slide5": 0.0, "slide5b2": 0.0, "slide5b3": 0.0, "slide6": 0.0}                               |
| tumour     | LCN2              | xgb/genes      |         0     | 0.000 [0.000, 0.000]    |  0     | -0.414 [-0.588, -0.291]  | {"slide1": 0.0, "slide2": 0.0, "slide3": 0.0, "slide4": 0.0, "slide5": 0.0, "slide5b2": 0.0, "slide5b3": 0.0, "slide6": 0.0}                               |
| tumour     | PIGR              | xgb/genes      |         0     | 0.000 [0.000, 0.000]    |  0     | -0.270 [-0.345, -0.160]  | {"slide1": 0.0, "slide2": 0.0, "slide3": 0.0, "slide4": 0.0, "slide5": 0.0, "slide5b2": 0.0, "slide5b3": 0.0, "slide6": 0.0}                               |
| tumour     | OLFM4             | xgb/genes      |         0     | 0.000 [0.000, 0.000]    |  0     | -0.244 [-0.291, -0.132]  | {"slide1": 0.0, "slide2": 0.0, "slide3": 0.0, "slide4": 0.0, "slide5": 0.0, "slide5b2": 0.0, "slide5b3": 0.0, "slide6": 0.0}                               |
| tumour     | SPINK1            | xgb/genes      |         0     | 0.000 [0.000, 0.000]    |  0     | -0.458 [-0.514, -0.331]  | {"slide1": 0.0, "slide2": 0.0, "slide3": 0.0, "slide4": 0.0, "slide5": 0.0, "slide5b2": 0.0, "slide5b3": 0.0, "slide6": 0.0}                               |
| tumour     | CEACAM1           | xgb/genes      |         0     | 0.000 [0.000, 0.000]    |  0     | 0.005 [-0.043, 0.052]    | {"slide1": 0.0, "slide2": 0.0, "slide3": 0.0, "slide4": 0.0, "slide5": 0.0, "slide5b2": 0.0, "slide5b3": 0.0, "slide6": 0.0}                               |
| tumour     | KRT6A             | xgb/genes      |         0     | 0.000 [0.000, 0.000]    |  0     | -0.314 [-0.379, -0.233]  | {"slide1": 0.0, "slide2": 0.0, "slide3": 0.0, "slide4": 0.0, "slide5": 0.0, "slide5b2": 0.0, "slide5b3": 0.0, "slide6": 0.0}                               |
| tumour     | KRT14             | xgb/genes      |         0     | 0.000 [0.000, 0.000]    |  0     | -0.275 [-0.342, -0.183]  | {"slide1": 0.0, "slide2": 0.0, "slide3": 0.0, "slide4": 0.0, "slide5": 0.0, "slide5b2": 0.0, "slide5b3": 0.0, "slide6": 0.0}                               |
| tumour     | KRT17             | xgb/genes      |         0     | 0.000 [0.000, 0.000]    |  0     | -0.471 [-0.507, -0.411]  | {"slide1": 0.0, "slide2": 0.0, "slide3": 0.0, "slide4": 0.0, "slide5": 0.0, "slide5b2": 0.0, "slide5b3": 0.0, "slide6": 0.0}                               |
| tumour     | S100A2            | xgb/genes      |         0     | 0.000 [0.000, 0.000]    |  0     | -0.464 [-0.592, -0.280]  | {"slide1": 0.0, "slide2": 0.0, "slide3": 0.0, "slide4": 0.0, "slide5": 0.0, "slide5b2": 0.0, "slide5b3": 0.0, "slide6": 0.0}                               |
| tumour     | SERPINB5          | xgb/genes      |         0     | 0.000 [0.000, 0.000]    |  0     | -0.183 [-0.236, -0.125]  | {"slide1": 0.0, "slide2": 0.0, "slide3": 0.0, "slide4": 0.0, "slide5": 0.0, "slide5b2": 0.0, "slide5b3": 0.0, "slide6": 0.0}                               |
| tumour     | AREG              | xgb/genes      |         0     | 0.000 [0.000, 0.000]    |  0     | -0.137 [-0.184, -0.101]  | {"slide1": 0.0, "slide2": 0.0, "slide3": 0.0, "slide4": 0.0, "slide5": 0.0, "slide5b2": 0.0, "slide5b3": 0.0, "slide6": 0.0}                               |
| tumour     | SLC2A1            | xgb/genes      |         0     | 0.000 [0.000, 0.000]    |  0     | -0.159 [-0.223, -0.120]  | {"slide1": 0.0, "slide2": 0.0, "slide3": 0.0, "slide4": 0.0, "slide5": 0.0, "slide5b2": 0.0, "slide5b3": 0.0, "slide6": 0.0}                               |
| tumour     | KRT15             | xgb/genes      |         0     | 0.000 [0.000, 0.000]    |  0     | 0.022 [-0.006, 0.057]    | {"slide1": 0.0, "slide2": 0.0, "slide3": 0.0, "slide4": 0.0, "slide5": 0.0, "slide5b2": 0.0, "slide5b3": 0.0, "slide6": 0.0}                               |
| tumour     | ITGA3             | xgb/genes      |         0     | 0.000 [0.000, 0.000]    |  0     | -0.074 [-0.167, 0.001]   | {"slide1": 0.0, "slide2": 0.0, "slide3": 0.0, "slide4": 0.0, "slide5": 0.0, "slide5b2": 0.0, "slide5b3": 0.0, "slide6": 0.0}                               |
| tumour     | MET               | xgb/genes      |         0     | 0.000 [0.000, 0.000]    |  0     | -0.117 [-0.154, -0.089]  | {"slide1": 0.0, "slide2": 0.0, "slide3": 0.0, "slide4": 0.0, "slide5": 0.0, "slide5b2": 0.0, "slide5b3": 0.0, "slide6": 0.0}                               |
| tumour     | genes_pooled      | xgb/genes      |         0     | 0.000 [0.000, 0.000]    |  0     | -0.286 [-0.327, -0.259]  | {"slide1": 0.0, "slide2": 0.0, "slide3": 0.0, "slide4": 0.0, "slide5": 0.0, "slide5b2": 0.0, "slide5b3": 0.0, "slide6": 0.0}                               |
| tumour     | tumour:classical  | ridge/own      |         0.618 | -0.074 [-0.346, 0.266]  |  0.692 |                          | {"slide1": 0.4658, "slide2": 0.2924, "slide3": -0.1959, "slide4": -0.3676, "slide5": -0.5256, "slide5b2": -0.5526, "slide5b3": -0.2652, "slide6": -0.2626} |
| tumour     | tumour:basal_like | ridge/own      |         0.397 | 0.375 [0.290, 0.433]    |  0.022 |                          | {"slide1": 0.4345, "slide2": 0.4959, "slide3": 0.3838, "slide4": 0.2253, "slide5": 0.4623, "slide5b2": 0.3251, "slide5b3": 0.3755, "slide6": 0.1833}       |
| tumour     | CEACAM6           | ridge/own      |         0.627 | 0.285 [0.192, 0.411]    |  0.342 |                          | {"slide1": 0.5384, "slide2": 0.1781, "slide3": 0.1063, "slide4": 0.3587, "slide5": 0.1601, "slide5b2": -0.014, "slide5b3": 0.3868, "slide6": 0.2954}       |
| tumour     | AGR2              | ridge/own      |         0.451 | 0.275 [0.224, 0.352]    |  0.176 |                          | {"slide1": 0.2815, "slide2": 0.332, "slide3": 0.3768, "slide4": 0.3562, "slide5": 0.5032, "slide5b2": 0.1271, "slide5b3": -0.0197, "slide6": 0.124}        |
| tumour     | S100P             | ridge/own      |         0.466 | 0.323 [0.250, 0.412]    |  0.143 |                          | {"slide1": 0.5415, "slide2": 0.2846, "slide3": 0.2382, "slide4": 0.2278, "slide5": 0.5017, "slide5b2": 0.3428, "slide5b3": 0.154, "slide6": 0.2749}        |
| tumour     | LCN2              | ridge/own      |         0.54  | 0.414 [0.291, 0.588]    |  0.126 |                          | {"slide1": 0.5202, "slide2": 0.4794, "slide3": 0.2864, "slide4": 0.7375, "slide5": 0.2035, "slide5b2": 0.2123, "slide5b3": 0.2799, "slide6": 0.4282}       |
| tumour     | PIGR              | ridge/own      |         0.535 | 0.270 [0.160, 0.345]    |  0.265 |                          | {"slide1": 0.2463, "slide2": 0.0902, "slide3": 0.4266, "slide4": 0.1372, "slide5": 0.4815, "slide5b2": 0.2936, "slide5b3": 0.0837, "slide6": 0.2629}       |
| tumour     | OLFM4             | ridge/own      |         0.604 | 0.244 [0.132, 0.291]    |  0.36  |                          | {"slide1": 0.1688, "slide2": 0.1406, "slide3": 0.0611, "slide4": 0.3084, "slide5": 0.4714, "slide5b2": 0.2915, "slide5b3": 0.1234, "slide6": 0.287}        |
| tumour     | SPINK1            | ridge/own      |         0.579 | 0.458 [0.331, 0.514]    |  0.121 |                          | {"slide1": 0.5589, "slide2": 0.153, "slide3": 0.5318, "slide4": 0.4459, "slide5": 0.6892, "slide5b2": 0.4064, "slide5b3": 0.3815, "slide6": 0.4855}        |
| tumour     | CEACAM1           | ridge/own      |         0.132 | -0.005 [-0.052, 0.043]  |  0.138 |                          | {"slide1": 0.0499, "slide2": -0.0006, "slide3": -0.0982, "slide4": -0.0611, "slide5": 0.0862, "slide5b2": -0.0091, "slide5b3": -0.1272, "slide6": 0.0868}  |
| tumour     | KRT6A             | ridge/own      |         0.574 | 0.314 [0.233, 0.379]    |  0.259 |                          | {"slide1": 0.3899, "slide2": 0.1483, "slide3": 0.425, "slide4": 0.2815, "slide5": 0.3937, "slide5b2": 0.2319, "slide5b3": 0.2031, "slide6": 0.2838}        |
| tumour     | KRT14             | ridge/own      |         0.307 | 0.275 [0.183, 0.342]    |  0.032 |                          | {"slide1": 0.1809, "slide2": 0.2067, "slide3": 0.1377, "slide4": 0.4004, "slide5": 0.4365, "slide5b2": 0.3205, "slide5b3": 0.3287, "slide6": 0.2539}       |
| tumour     | KRT17             | ridge/own      |         0.428 | 0.471 [0.411, 0.507]    | -0.043 |                          | {"slide1": 0.5442, "slide2": 0.5138, "slide3": 0.2532, "slide4": 0.4454, "slide5": 0.4415, "slide5b2": 0.3098, "slide5b3": 0.6051, "slide6": 0.4446}       |
| tumour     | S100A2            | ridge/own      |         0.453 | 0.464 [0.280, 0.592]    | -0.011 |                          | {"slide1": 0.6555, "slide2": 0.3394, "slide3": 0.2367, "slide4": 0.4076, "slide5": 0.3436, "slide5b2": 0.3385, "slide5b3": 0.4064, "slide6": 0.2432}       |
| tumour     | SERPINB5          | ridge/own      |         0.147 | 0.183 [0.125, 0.236]    | -0.036 |                          | {"slide1": 0.1647, "slide2": 0.2514, "slide3": 0.0534, "slide4": 0.1351, "slide5": 0.113, "slide5b2": 0.2736, "slide5b3": 0.4468, "slide6": 0.1623}        |
| tumour     | AREG              | ridge/own      |         0.157 | 0.137 [0.101, 0.184]    |  0.021 |                          | {"slide1": 0.1672, "slide2": 0.256, "slide3": 0.1824, "slide4": 0.1155, "slide5": 0.0957, "slide5b2": 0.0807, "slide5b3": 0.0882, "slide6": 0.0906}        |
| tumour     | SLC2A1            | ridge/own      |         0.217 | 0.159 [0.120, 0.223]    |  0.059 |                          | {"slide1": 0.0916, "slide2": 0.3208, "slide3": 0.1968, "slide4": 0.1015, "slide5": 0.1737, "slide5b2": 0.0651, "slide5b3": 0.2005, "slide6": 0.1612}       |
| tumour     | KRT15             | ridge/own      |         0.089 | -0.022 [-0.057, 0.006]  |  0.112 |                          | {"slide1": -0.0725, "slide2": -0.0247, "slide3": -0.0189, "slide4": -0.0607, "slide5": -0.0565, "slide5b2": 0.0665, "slide5b3": 0.0528, "slide6": -0.0464} |
| tumour     | ITGA3             | ridge/own      |         0.238 | 0.074 [-0.001, 0.167]   |  0.164 |                          | {"slide1": 0.1895, "slide2": 0.2637, "slide3": -0.1041, "slide4": 0.0313, "slide5": 0.0184, "slide5b2": -0.0186, "slide5b3": 0.1241, "slide6": -0.0019}    |
| tumour     | MET               | ridge/own      |         0.086 | 0.117 [0.089, 0.154]    | -0.03  |                          | {"slide1": 0.1229, "slide2": 0.1679, "slide3": 0.0458, "slide4": 0.1338, "slide5": 0.0486, "slide5b2": 0.0699, "slide5b3": 0.2043, "slide6": 0.1632}       |
| tumour     | genes_pooled      | ridge/own      |         0.425 | 0.286 [0.259, 0.327]    |  0.139 |                          | {"slide1": 0.3834, "slide2": 0.259, "slide3": 0.2377, "slide4": 0.3025, "slide5": 0.3637, "slide5b2": 0.2107, "slide5b3": 0.2638, "slide6": 0.2585}        |
| tumour     | tumour:classical  | ridge/programs |         0.693 | 0.471 [0.356, 0.616]    |  0.221 | 0.546 [0.357, 0.708]     | {"slide1": 0.7376, "slide2": 0.5729, "slide3": 0.4411, "slide4": 0.3048, "slide5": 0.3209, "slide5b2": 0.328, "slide5b3": 0.3836, "slide6": 0.3508}        |
| tumour     | tumour:basal_like | ridge/programs |         0.438 | 0.423 [0.348, 0.469]    |  0.015 | 0.047 [0.027, 0.070]     | {"slide1": 0.4738, "slide2": 0.4995, "slide3": 0.4761, "slide4": 0.3021, "slide5": 0.5111, "slide5b2": 0.3772, "slide5b3": 0.3997, "slide6": 0.2294}       |
| tumour     | CEACAM6           | ridge/programs |         0.654 | 0.509 [0.418, 0.600]    |  0.145 | 0.224 [0.159, 0.280]     | {"slide1": 0.6658, "slide2": 0.3083, "slide3": 0.4733, "slide4": 0.5756, "slide5": 0.4698, "slide5b2": 0.3414, "slide5b3": 0.5685, "slide6": 0.5536}       |
| tumour     | AGR2              | ridge/programs |         0.461 | 0.400 [0.331, 0.459]    |  0.061 | 0.125 [0.082, 0.148]     | {"slide1": 0.3862, "slide2": 0.388, "slide3": 0.5151, "slide4": 0.4514, "slide5": 0.6054, "slide5b2": 0.3088, "slide5b3": 0.192, "slide6": 0.211}          |
| tumour     | S100P             | ridge/programs |         0.485 | 0.440 [0.364, 0.503]    |  0.045 | 0.118 [0.083, 0.139]     | {"slide1": 0.6035, "slide2": 0.3548, "slide3": 0.3809, "slide4": 0.3507, "slide5": 0.5841, "slide5b2": 0.496, "slide5b3": 0.3308, "slide6": 0.4058}        |
| tumour     | LCN2              | ridge/programs |         0.556 | 0.472 [0.366, 0.624]    |  0.084 | 0.058 [0.030, 0.077]     | {"slide1": 0.5793, "slide2": 0.5006, "slide3": 0.3517, "slide4": 0.7518, "slide5": 0.2636, "slide5b2": 0.2991, "slide5b3": 0.3734, "slide6": 0.5013}       |
| tumour     | PIGR              | ridge/programs |         0.549 | 0.373 [0.266, 0.457]    |  0.177 | 0.103 [0.093, 0.119]     | {"slide1": 0.382, "slide2": 0.1901, "slide3": 0.5407, "slide4": 0.2369, "slide5": 0.5499, "slide5b2": 0.3938, "slide5b3": 0.1865, "slide6": 0.3798}        |
| tumour     | OLFM4             | ridge/programs |         0.616 | 0.390 [0.285, 0.437]    |  0.226 | 0.146 [0.133, 0.177]     | {"slide1": 0.3618, "slide2": 0.2668, "slide3": 0.2332, "slide4": 0.4261, "slide5": 0.5576, "slide5b2": 0.4367, "slide5b3": 0.299, "slide6": 0.4729}        |
| tumour     | SPINK1            | ridge/programs |         0.593 | 0.552 [0.416, 0.615]    |  0.041 | 0.094 [0.082, 0.107]     | {"slide1": 0.6234, "slide2": 0.2217, "slide3": 0.6463, "slide4": 0.5466, "slide5": 0.7493, "slide5b2": 0.514, "slide5b3": 0.488, "slide6": 0.6012}         |
| tumour     | CEACAM1           | ridge/programs |         0.138 | 0.090 [0.035, 0.139]    |  0.048 | 0.095 [0.069, 0.112]     | {"slide1": 0.1201, "slide2": 0.0539, "slide3": 0.0301, "slide4": 0.0035, "slide5": 0.147, "slide5b2": 0.1014, "slide5b3": 0.06, "slide6": 0.1874}          |
| tumour     | KRT6A             | ridge/programs |         0.588 | 0.371 [0.286, 0.442]    |  0.217 | 0.057 [0.044, 0.070]     | {"slide1": 0.4317, "slide2": 0.211, "slide3": 0.5105, "slide4": 0.3116, "slide5": 0.4239, "slide5b2": 0.2932, "slide5b3": 0.2829, "slide6": 0.341}         |
| tumour     | KRT14             | ridge/programs |         0.301 | 0.264 [0.175, 0.330]    |  0.037 | -0.011 [-0.019, 0.001]   | {"slide1": 0.179, "slide2": 0.1857, "slide3": 0.1184, "slide4": 0.4031, "slide5": 0.4727, "slide5b2": 0.287, "slide5b3": 0.2696, "slide6": 0.2663}         |
| tumour     | KRT17             | ridge/programs |         0.433 | 0.504 [0.455, 0.534]    | -0.071 | 0.033 [0.020, 0.046]     | {"slide1": 0.5743, "slide2": 0.5217, "slide3": 0.325, "slide4": 0.4916, "slide5": 0.4736, "slide5b2": 0.3693, "slide5b3": 0.626, "slide6": 0.4768}         |
| tumour     | S100A2            | ridge/programs |         0.456 | 0.473 [0.307, 0.599]    | -0.016 | 0.009 [-0.001, 0.030]    | {"slide1": 0.6584, "slide2": 0.3427, "slide3": 0.2995, "slide4": 0.4124, "slide5": 0.3498, "slide5b2": 0.3355, "slide5b3": 0.3907, "slide6": 0.2501}       |
| tumour     | SERPINB5          | ridge/programs |         0.153 | 0.143 [0.079, 0.188]    |  0.009 | -0.040 [-0.068, -0.009]  | {"slide1": 0.1556, "slide2": 0.2366, "slide3": -0.0353, "slide4": 0.1048, "slide5": 0.0068, "slide5b2": 0.2006, "slide5b3": 0.393, "slide6": 0.1643}       |
| tumour     | AREG              | ridge/programs |         0.149 | 0.138 [0.101, 0.185]    |  0.011 | 0.002 [-0.001, 0.003]    | {"slide1": 0.1716, "slide2": 0.253, "slide3": 0.1868, "slide4": 0.1151, "slide5": 0.1011, "slide5b2": 0.0797, "slide5b3": 0.0857, "slide6": 0.0881}        |
| tumour     | SLC2A1            | ridge/programs |         0.224 | 0.198 [0.148, 0.246]    |  0.026 | 0.039 [0.012, 0.063]     | {"slide1": 0.1228, "slide2": 0.3264, "slide3": 0.2112, "slide4": 0.1285, "slide5": 0.2216, "slide5b2": 0.1387, "slide5b3": 0.3156, "slide6": 0.1701}       |
| tumour     | KRT15             | ridge/programs |         0.079 | -0.003 [-0.046, 0.027]  |  0.082 | 0.019 [0.001, 0.045]     | {"slide1": -0.089, "slide2": -0.0297, "slide3": 0.0675, "slide4": -0.0271, "slide5": -0.0008, "slide5b2": 0.0663, "slide5b3": 0.0243, "slide6": -0.0014}   |
| tumour     | ITGA3             | ridge/programs |         0.239 | 0.194 [0.129, 0.259]    |  0.045 | 0.120 [0.088, 0.141]     | {"slide1": 0.2818, "slide2": 0.3299, "slide3": 0.0491, "slide4": 0.143, "slide5": 0.1567, "slide5b2": 0.1374, "slide5b3": 0.2664, "slide6": 0.1133}        |
| tumour     | MET               | ridge/programs |         0.09  | 0.131 [0.109, 0.162]    | -0.041 | 0.014 [-0.003, 0.028]    | {"slide1": 0.1056, "slide2": 0.1459, "slide3": 0.0932, "slide4": 0.1455, "slide5": 0.0781, "slide5b2": 0.0882, "slide5b3": 0.2008, "slide6": 0.1872}       |
| tumour     | genes_pooled      | ridge/programs |         0.435 | 0.372 [0.339, 0.399]    |  0.063 | 0.086 [0.065, 0.100]     | {"slide1": 0.4418, "slide2": 0.3102, "slide3": 0.3531, "slide4": 0.3786, "slide5": 0.442, "slide5b2": 0.3243, "slide5b3": 0.3673, "slide6": 0.3531}        |
| tumour     | tumour:classical  | graph/sage     |         0.597 | 0.466 [0.351, 0.593]    |  0.131 | 0.540 [0.325, 0.724]     | {"slide1": 0.7125, "slide2": 0.536, "slide3": 0.4297, "slide4": 0.2038, "slide5": 0.2368, "slide5b2": 0.4222, "slide5b3": 0.4489, "slide6": 0.3716}        |
| tumour     | tumour:basal_like | graph/sage     |         0.281 | 0.449 [0.358, 0.504]    | -0.168 | 0.073 [0.024, 0.106]     | {"slide1": 0.4167, "slide2": 0.5552, "slide3": 0.4689, "slide4": 0.263, "slide5": 0.537, "slide5b2": 0.425, "slide5b3": 0.5469, "slide6": 0.3194}          |
| tumour     | CEACAM6           | graph/sage     |         0.651 | 0.475 [0.374, 0.562]    |  0.176 | 0.191 [0.133, 0.240]     | {"slide1": 0.6715, "slide2": 0.3775, "slide3": 0.281, "slide4": 0.4899, "slide5": 0.4157, "slide5b2": 0.2992, "slide5b3": 0.6284, "slide6": 0.4119}        |
| tumour     | AGR2              | graph/sage     |         0.465 | 0.329 [0.281, 0.414]    |  0.136 | 0.054 [-0.009, 0.115]    | {"slide1": 0.4831, "slide2": 0.4847, "slide3": 0.3353, "slide4": 0.3067, "slide5": 0.4437, "slide5b2": 0.1818, "slide5b3": 0.2229, "slide6": 0.1783}       |
| tumour     | S100P             | graph/sage     |         0.474 | 0.419 [0.292, 0.514]    |  0.055 | 0.096 [0.027, 0.129]     | {"slide1": 0.6528, "slide2": 0.3784, "slide3": 0.3745, "slide4": 0.1731, "slide5": 0.5749, "slide5b2": 0.4783, "slide5b3": 0.3389, "slide6": 0.323}        |
| tumour     | LCN2              | graph/sage     |         0.554 | 0.477 [0.379, 0.611]    |  0.077 | 0.063 [0.012, 0.094]     | {"slide1": 0.6112, "slide2": 0.4848, "slide3": 0.3778, "slide4": 0.7228, "slide5": 0.2909, "slide5b2": 0.3123, "slide5b3": 0.3944, "slide6": 0.4775}       |
| tumour     | PIGR              | graph/sage     |         0.554 | 0.303 [0.207, 0.340]    |  0.251 | 0.033 [-0.028, 0.080]    | {"slide1": 0.3671, "slide2": 0.2007, "slide3": 0.3555, "slide4": 0.1793, "slide5": 0.5446, "slide5b2": 0.3873, "slide5b3": 0.113, "slide6": 0.2109}        |
| tumour     | OLFM4             | graph/sage     |         0.648 | 0.483 [0.374, 0.533]    |  0.165 | 0.240 [0.213, 0.290]     | {"slide1": 0.4781, "slide2": 0.32, "slide3": 0.357, "slide4": 0.5213, "slide5": 0.6441, "slide5b2": 0.5144, "slide5b3": 0.4069, "slide6": 0.5831}          |
| tumour     | SPINK1            | graph/sage     |         0.603 | 0.469 [0.375, 0.549]    |  0.133 | 0.011 [-0.017, 0.073]    | {"slide1": 0.7077, "slide2": 0.2309, "slide3": 0.4937, "slide4": 0.4593, "slide5": 0.6452, "slide5b2": 0.3942, "slide5b3": 0.3712, "slide6": 0.5122}       |
| tumour     | CEACAM1           | graph/sage     |         0.139 | 0.177 [0.100, 0.243]    | -0.038 | 0.183 [0.112, 0.240]     | {"slide1": 0.3981, "slide2": 0.1795, "slide3": 0.0835, "slide4": 0.0606, "slide5": 0.1694, "slide5b2": 0.193, "slide5b3": 0.3575, "slide6": 0.1427}        |
| tumour     | KRT6A             | graph/sage     |         0.59  | 0.379 [0.302, 0.431]    |  0.211 | 0.065 [0.033, 0.090]     | {"slide1": 0.4484, "slide2": 0.2526, "slide3": 0.4664, "slide4": 0.2654, "slide5": 0.5017, "slide5b2": 0.3057, "slide5b3": 0.3152, "slide6": 0.3633}       |
| tumour     | KRT14             | graph/sage     |         0.307 | 0.398 [0.279, 0.489]    | -0.091 | 0.123 [0.084, 0.155]     | {"slide1": 0.2502, "slide2": 0.2907, "slide3": 0.2793, "slide4": 0.5061, "slide5": 0.537, "slide5b2": 0.4801, "slide5b3": 0.5694, "slide6": 0.3663}        |
| tumour     | KRT17             | graph/sage     |         0.413 | 0.485 [0.440, 0.529]    | -0.072 | 0.014 [-0.000, 0.045]    | {"slide1": 0.5879, "slide2": 0.5162, "slide3": 0.36, "slide4": 0.4466, "slide5": 0.4648, "slide5b2": 0.3639, "slide5b3": 0.5581, "slide6": 0.4448}         |
| tumour     | S100A2            | graph/sage     |         0.442 | 0.465 [0.341, 0.546]    | -0.023 | 0.001 [-0.045, 0.063]    | {"slide1": 0.5885, "slide2": 0.388, "slide3": 0.3253, "slide4": 0.4294, "slide5": 0.3881, "slide5b2": 0.369, "slide5b3": 0.4619, "slide6": 0.2835}         |
| tumour     | SERPINB5          | graph/sage     |         0.148 | 0.200 [0.148, 0.248]    | -0.052 | 0.017 [0.001, 0.042]     | {"slide1": 0.1933, "slide2": 0.3226, "slide3": 0.0963, "slide4": 0.1468, "slide5": 0.1695, "slide5b2": 0.2361, "slide5b3": 0.3802, "slide6": 0.169}        |
| tumour     | AREG              | graph/sage     |         0.146 | 0.141 [0.114, 0.198]    |  0.006 | 0.004 [-0.027, 0.034]    | {"slide1": 0.1796, "slide2": 0.3403, "slide3": 0.1264, "slide4": 0.1245, "slide5": 0.0865, "slide5b2": 0.0869, "slide5b3": 0.1661, "slide6": 0.1016}       |
| tumour     | SLC2A1            | graph/sage     |         0.236 | 0.320 [0.231, 0.393]    | -0.083 | 0.161 [0.090, 0.206]     | {"slide1": 0.3166, "slide2": 0.5058, "slide3": 0.3725, "slide4": 0.1795, "slide5": 0.2652, "slide5b2": 0.2428, "slide5b3": 0.5712, "slide6": 0.1907}       |
| tumour     | KRT15             | graph/sage     |         0.083 | 0.145 [0.068, 0.205]    | -0.061 | 0.167 [0.111, 0.203]     | {"slide1": 0.1272, "slide2": 0.0896, "slide3": 0.1594, "slide4": 0.0523, "slide5": 0.0537, "slide5b2": 0.2615, "slide5b3": 0.4035, "slide6": 0.0427}       |
| tumour     | ITGA3             | graph/sage     |         0.227 | 0.352 [0.263, 0.413]    | -0.125 | 0.278 [0.207, 0.323]     | {"slide1": 0.4715, "slide2": 0.459, "slide3": 0.287, "slide4": 0.2512, "slide5": 0.1972, "slide5b2": 0.2783, "slide5b3": 0.5748, "slide6": 0.1774}         |
| tumour     | MET               | graph/sage     |         0.076 | 0.174 [0.150, 0.210]    | -0.097 | 0.057 [0.037, 0.076]     | {"slide1": 0.2199, "slide2": 0.2507, "slide3": 0.1123, "slide4": 0.1724, "slide5": 0.0751, "slide5b2": 0.1353, "slide5b3": 0.3101, "slide6": 0.181}        |
| tumour     | genes_pooled      | graph/sage     |         0.436 | 0.385 [0.336, 0.427]    |  0.051 | 0.099 [0.070, 0.115]     | {"slide1": 0.4965, "slide2": 0.3637, "slide3": 0.3254, "slide4": 0.3495, "slide5": 0.433, "slide5b2": 0.3395, "slide5b3": 0.4195, "slide6": 0.3246}        |

## Phase 5 permutation importance (top groups)

| receiver   | target       | group                         |   r2_drop |   ci_low |   ci_high |   n_tumours |
|:-----------|:-------------|:------------------------------|----------:|---------:|----------:|------------:|
| myeloid    | genes_pooled | genes:other_G_genes           |    0.138  |   0.1074 |    0.1758 |          13 |
| myeloid    | genes_pooled | class:myeloid                 |    0.0332 |   0.022  |    0.0412 |          13 |
| myeloid    | genes_pooled | genes:myeloid:SPP1_TAM        |    0.0325 |   0.0164 |    0.0474 |          13 |
| myeloid    | genes_pooled | genes:myeloid:C1Q_TAM         |    0.0242 |   0.0111 |    0.0398 |          13 |
| myeloid    | genes_pooled | class:CAF_fibroblast          |    0.0107 |   0.0092 |    0.0123 |          13 |
| myeloid    | genes_pooled | class:epithelial_ductal_tumor |    0.0041 |   0.0013 |    0.0083 |          13 |

## Phase 5 in-silico perturbation (model-based hypotheses)

| receiver   | target       | neighbour_program   |   quantile |   predicted_change |   ci_low |   ci_high | label                  |
|:-----------|:-------------|:--------------------|-----------:|-------------------:|---------:|----------:|:-----------------------|
| myeloid    | genes_pooled | tumour:classical    |        0.1 |            -0.0025 |  -0.0068 |    0.0005 | model-based hypothesis |
| myeloid    | C1QA         | tumour:classical    |        0.1 |            -0.0004 |  -0.006  |    0.0041 | model-based hypothesis |
| myeloid    | C1QB         | tumour:classical    |        0.1 |            -0.0072 |  -0.0158 |   -0.0014 | model-based hypothesis |
| myeloid    | C1QC         | tumour:classical    |        0.1 |            -0.0057 |  -0.0147 |    0.0004 | model-based hypothesis |
| myeloid    | CD163        | tumour:classical    |        0.1 |             0.005  |   0.003  |    0.0071 | model-based hypothesis |
| myeloid    | MRC1         | tumour:classical    |        0.1 |             0.0057 |   0.0032 |    0.0089 | model-based hypothesis |
| myeloid    | SELENOP      | tumour:classical    |        0.1 |             0.002  |  -0.0002 |    0.005  | model-based hypothesis |
| myeloid    | SPP1         | tumour:classical    |        0.1 |            -0.0195 |  -0.0348 |   -0.0068 | model-based hypothesis |
| myeloid    | GPNMB        | tumour:classical    |        0.1 |            -0.0008 |  -0.0057 |    0.0028 | model-based hypothesis |
| myeloid    | CD9          | tumour:classical    |        0.1 |            -0.0032 |  -0.0073 |   -0.0001 | model-based hypothesis |
| myeloid    | FABP5        | tumour:classical    |        0.1 |            -0.0001 |  -0.0025 |    0.0019 | model-based hypothesis |
| myeloid    | LGALS3       | tumour:classical    |        0.1 |            -0.0033 |  -0.0082 |   -0      | model-based hypothesis |
| myeloid    | genes_pooled | tumour:classical    |        0.9 |            -0.0038 |  -0.007  |   -0.0013 | model-based hypothesis |
| myeloid    | C1QA         | tumour:classical    |        0.9 |            -0.0059 |  -0.0119 |   -0.0013 | model-based hypothesis |
| myeloid    | C1QB         | tumour:classical    |        0.9 |             0.0028 |   0.0009 |    0.0047 | model-based hypothesis |
| myeloid    | C1QC         | tumour:classical    |        0.9 |             0.0018 |  -0.002  |    0.005  | model-based hypothesis |
| myeloid    | CD163        | tumour:classical    |        0.9 |            -0.0112 |  -0.0153 |   -0.0076 | model-based hypothesis |
| myeloid    | MRC1         | tumour:classical    |        0.9 |            -0.0138 |  -0.0236 |   -0.0066 | model-based hypothesis |
| myeloid    | SELENOP      | tumour:classical    |        0.9 |            -0.0055 |  -0.0119 |   -0.0002 | model-based hypothesis |
| myeloid    | SPP1         | tumour:classical    |        0.9 |             0.0044 |  -0.0073 |    0.0165 | model-based hypothesis |
| myeloid    | GPNMB        | tumour:classical    |        0.9 |            -0.0055 |  -0.0095 |   -0.0022 | model-based hypothesis |
| myeloid    | CD9          | tumour:classical    |        0.9 |            -0.001  |  -0.0035 |    0.0014 | model-based hypothesis |
| myeloid    | FABP5        | tumour:classical    |        0.9 |            -0.0078 |  -0.0134 |   -0.0034 | model-based hypothesis |
| myeloid    | LGALS3       | tumour:classical    |        0.9 |            -0.0003 |  -0.0036 |    0.0024 | model-based hypothesis |
| myeloid    | genes_pooled | tumour:basal_like   |        0.1 |            -0.0083 |  -0.0122 |   -0.0052 | model-based hypothesis |
| myeloid    | C1QA         | tumour:basal_like   |        0.1 |            -0.007  |  -0.0111 |   -0.0034 | model-based hypothesis |
| myeloid    | C1QB         | tumour:basal_like   |        0.1 |            -0.0005 |  -0.0074 |    0.0085 | model-based hypothesis |
| myeloid    | C1QC         | tumour:basal_like   |        0.1 |            -0.0064 |  -0.0119 |   -0.0007 | model-based hypothesis |
| myeloid    | CD163        | tumour:basal_like   |        0.1 |             0.0009 |  -0.0026 |    0.0053 | model-based hypothesis |
| myeloid    | MRC1         | tumour:basal_like   |        0.1 |            -0.005  |  -0.0103 |   -0.0016 | model-based hypothesis |
| myeloid    | SELENOP      | tumour:basal_like   |        0.1 |             0.0088 |   0.0027 |    0.0183 | model-based hypothesis |
| myeloid    | SPP1         | tumour:basal_like   |        0.1 |            -0.0473 |  -0.0845 |   -0.0194 | model-based hypothesis |
| myeloid    | GPNMB        | tumour:basal_like   |        0.1 |            -0.0138 |  -0.0209 |   -0.008  | model-based hypothesis |
| myeloid    | CD9          | tumour:basal_like   |        0.1 |            -0.009  |  -0.0181 |   -0.0029 | model-based hypothesis |
| myeloid    | FABP5        | tumour:basal_like   |        0.1 |            -0.0057 |  -0.0113 |   -0.0021 | model-based hypothesis |
| myeloid    | LGALS3       | tumour:basal_like   |        0.1 |            -0.0066 |  -0.0108 |   -0.0034 | model-based hypothesis |
| myeloid    | genes_pooled | tumour:basal_like   |        0.9 |             0.008  |   0.0041 |    0.013  | model-based hypothesis |
| myeloid    | C1QA         | tumour:basal_like   |        0.9 |             0.0082 |   0.0018 |    0.0159 | model-based hypothesis |
| myeloid    | C1QB         | tumour:basal_like   |        0.9 |             0.0018 |  -0.0096 |    0.0119 | model-based hypothesis |
| myeloid    | C1QC         | tumour:basal_like   |        0.9 |             0.0061 |  -0.002  |    0.015  | model-based hypothesis |
| myeloid    | CD163        | tumour:basal_like   |        0.9 |            -0.0016 |  -0.0077 |    0.0038 | model-based hypothesis |
| myeloid    | MRC1         | tumour:basal_like   |        0.9 |             0.0062 |   0.0023 |    0.0106 | model-based hypothesis |
| myeloid    | SELENOP      | tumour:basal_like   |        0.9 |            -0.0088 |  -0.0167 |   -0.0033 | model-based hypothesis |
| myeloid    | SPP1         | tumour:basal_like   |        0.9 |             0.0479 |   0.0195 |    0.0917 | model-based hypothesis |
| myeloid    | GPNMB        | tumour:basal_like   |        0.9 |             0.0116 |   0.0038 |    0.022  | model-based hypothesis |
| myeloid    | CD9          | tumour:basal_like   |        0.9 |             0.0045 |  -0.0016 |    0.0134 | model-based hypothesis |
| myeloid    | FABP5        | tumour:basal_like   |        0.9 |             0.0061 |   0.0032 |    0.0107 | model-based hypothesis |
| myeloid    | LGALS3       | tumour:basal_like   |        0.9 |             0.0061 |   0.0027 |    0.0115 | model-based hypothesis |
| myeloid    | genes_pooled | CAF:myCAF           |        0.1 |            -0.0097 |  -0.0141 |   -0.0057 | model-based hypothesis |
| myeloid    | C1QA         | CAF:myCAF           |        0.1 |            -0.0108 |  -0.0205 |   -0.0015 | model-based hypothesis |
| myeloid    | C1QB         | CAF:myCAF           |        0.1 |            -0.0121 |  -0.0228 |   -0.002  | model-based hypothesis |
| myeloid    | C1QC         | CAF:myCAF           |        0.1 |            -0.0252 |  -0.0354 |   -0.0158 | model-based hypothesis |
| myeloid    | CD163        | CAF:myCAF           |        0.1 |            -0.0278 |  -0.0429 |   -0.0137 | model-based hypothesis |
| myeloid    | MRC1         | CAF:myCAF           |        0.1 |            -0.0171 |  -0.0222 |   -0.0129 | model-based hypothesis |
| myeloid    | SELENOP      | CAF:myCAF           |        0.1 |            -0.0107 |  -0.0178 |   -0.0035 | model-based hypothesis |
| myeloid    | SPP1         | CAF:myCAF           |        0.1 |             0.008  |  -0.0011 |    0.0165 | model-based hypothesis |
| myeloid    | GPNMB        | CAF:myCAF           |        0.1 |            -0.0231 |  -0.0324 |   -0.0158 | model-based hypothesis |
| myeloid    | CD9          | CAF:myCAF           |        0.1 |            -0.0015 |  -0.004  |    0.0011 | model-based hypothesis |
| myeloid    | FABP5        | CAF:myCAF           |        0.1 |             0.0146 |   0.0107 |    0.0185 | model-based hypothesis |
| myeloid    | LGALS3       | CAF:myCAF           |        0.1 |            -0.0006 |  -0.0033 |    0.0019 | model-based hypothesis |
| myeloid    | genes_pooled | CAF:myCAF           |        0.9 |             0.0259 |   0.0206 |    0.0311 | model-based hypothesis |
| myeloid    | C1QA         | CAF:myCAF           |        0.9 |             0.0362 |   0.027  |    0.0459 | model-based hypothesis |
| myeloid    | C1QB         | CAF:myCAF           |        0.9 |             0.0275 |   0.0194 |    0.0365 | model-based hypothesis |
| myeloid    | C1QC         | CAF:myCAF           |        0.9 |             0.0628 |   0.0508 |    0.0754 | model-based hypothesis |
| myeloid    | CD163        | CAF:myCAF           |        0.9 |             0.0206 |   0.0048 |    0.039  | model-based hypothesis |
| myeloid    | MRC1         | CAF:myCAF           |        0.9 |             0.0221 |   0.0144 |    0.0305 | model-based hypothesis |
| myeloid    | SELENOP      | CAF:myCAF           |        0.9 |             0.0107 |   0.0004 |    0.0223 | model-based hypothesis |
| myeloid    | SPP1         | CAF:myCAF           |        0.9 |             0.0373 |   0.0147 |    0.0584 | model-based hypothesis |
| myeloid    | GPNMB        | CAF:myCAF           |        0.9 |             0.0677 |   0.0519 |    0.086  | model-based hypothesis |
| myeloid    | CD9          | CAF:myCAF           |        0.9 |             0.013  |   0.0066 |    0.0193 | model-based hypothesis |
| myeloid    | FABP5        | CAF:myCAF           |        0.9 |            -0.0241 |  -0.0303 |   -0.0186 | model-based hypothesis |
| myeloid    | LGALS3       | CAF:myCAF           |        0.9 |             0.0106 |   0.0044 |    0.0172 | model-based hypothesis |
| myeloid    | genes_pooled | CAF:CXCL12          |        0.1 |            -0.0034 |  -0.0054 |   -0.0018 | model-based hypothesis |
| myeloid    | C1QA         | CAF:CXCL12          |        0.1 |            -0.0062 |  -0.0096 |   -0.0035 | model-based hypothesis |
| myeloid    | C1QB         | CAF:CXCL12          |        0.1 |            -0.0032 |  -0.0059 |   -0.0008 | model-based hypothesis |
| myeloid    | C1QC         | CAF:CXCL12          |        0.1 |            -0.0055 |  -0.009  |   -0.0027 | model-based hypothesis |
| myeloid    | CD163        | CAF:CXCL12          |        0.1 |            -0.0111 |  -0.0163 |   -0.007  | model-based hypothesis |
| myeloid    | MRC1         | CAF:CXCL12          |        0.1 |            -0.0055 |  -0.0087 |   -0.0032 | model-based hypothesis |
| myeloid    | SELENOP      | CAF:CXCL12          |        0.1 |            -0.0198 |  -0.0276 |   -0.0132 | model-based hypothesis |
| myeloid    | SPP1         | CAF:CXCL12          |        0.1 |             0.0095 |   0.0058 |    0.0139 | model-based hypothesis |
| myeloid    | GPNMB        | CAF:CXCL12          |        0.1 |            -0.0007 |  -0.0036 |    0.0016 | model-based hypothesis |
| myeloid    | CD9          | CAF:CXCL12          |        0.1 |             0.0029 |   0.0019 |    0.0042 | model-based hypothesis |
| myeloid    | FABP5        | CAF:CXCL12          |        0.1 |             0.0023 |   0.0012 |    0.0036 | model-based hypothesis |
| myeloid    | LGALS3       | CAF:CXCL12          |        0.1 |             0.0003 |   0.0001 |    0.0005 | model-based hypothesis |
| myeloid    | genes_pooled | CAF:CXCL12          |        0.9 |             0.0128 |   0.0074 |    0.0198 | model-based hypothesis |
| myeloid    | C1QA         | CAF:CXCL12          |        0.9 |             0.0239 |   0.0149 |    0.0349 | model-based hypothesis |
| myeloid    | C1QB         | CAF:CXCL12          |        0.9 |             0.0133 |   0.0039 |    0.0236 | model-based hypothesis |
| myeloid    | C1QC         | CAF:CXCL12          |        0.9 |             0.0233 |   0.0124 |    0.0359 | model-based hypothesis |
| myeloid    | CD163        | CAF:CXCL12          |        0.9 |             0.0441 |   0.0312 |    0.0601 | model-based hypothesis |
| myeloid    | MRC1         | CAF:CXCL12          |        0.9 |             0.0222 |   0.014  |    0.0327 | model-based hypothesis |
| myeloid    | SELENOP      | CAF:CXCL12          |        0.9 |             0.074  |   0.0581 |    0.0937 | model-based hypothesis |
| myeloid    | SPP1         | CAF:CXCL12          |        0.9 |            -0.0412 |  -0.0519 |   -0.0329 | model-based hypothesis |
| myeloid    | GPNMB        | CAF:CXCL12          |        0.9 |             0.0007 |  -0.0064 |    0.0085 | model-based hypothesis |
| myeloid    | CD9          | CAF:CXCL12          |        0.9 |            -0.0104 |  -0.0124 |   -0.0084 | model-based hypothesis |
| myeloid    | FABP5        | CAF:CXCL12          |        0.9 |            -0.0075 |  -0.0088 |   -0.0061 | model-based hypothesis |
| myeloid    | LGALS3       | CAF:CXCL12          |        0.9 |            -0.0015 |  -0.0028 |   -0.0001 | model-based hypothesis |
| myeloid    | genes_pooled | T:cytotoxic         |        0.1 |            -0.0011 |  -0.0019 |   -0.0004 | model-based hypothesis |
| myeloid    | C1QA         | T:cytotoxic         |        0.1 |            -0.0038 |  -0.006  |   -0.0019 | model-based hypothesis |
| myeloid    | C1QB         | T:cytotoxic         |        0.1 |            -0.0026 |  -0.0041 |   -0.0013 | model-based hypothesis |
| myeloid    | C1QC         | T:cytotoxic         |        0.1 |            -0.0036 |  -0.0055 |   -0.0019 | model-based hypothesis |
| myeloid    | CD163        | T:cytotoxic         |        0.1 |            -0.0015 |  -0.0035 |    0.0003 | model-based hypothesis |
| myeloid    | MRC1         | T:cytotoxic         |        0.1 |            -0.0018 |  -0.0031 |   -0.0007 | model-based hypothesis |
| myeloid    | SELENOP      | T:cytotoxic         |        0.1 |            -0.0014 |  -0.0023 |   -0.0006 | model-based hypothesis |
| myeloid    | SPP1         | T:cytotoxic         |        0.1 |             0.0002 |  -0.0005 |    0.0011 | model-based hypothesis |
| myeloid    | GPNMB        | T:cytotoxic         |        0.1 |            -0.0001 |  -0.0009 |    0.0008 | model-based hypothesis |
| myeloid    | CD9          | T:cytotoxic         |        0.1 |             0.0009 |   0.0004 |    0.0017 | model-based hypothesis |
| myeloid    | FABP5        | T:cytotoxic         |        0.1 |             0.0011 |   0.0003 |    0.0022 | model-based hypothesis |
| myeloid    | LGALS3       | T:cytotoxic         |        0.1 |             0.0003 |  -0.0002 |    0.001  | model-based hypothesis |
| myeloid    | genes_pooled | T:cytotoxic         |        0.9 |             0.0042 |   0.0028 |    0.0057 | model-based hypothesis |
| myeloid    | C1QA         | T:cytotoxic         |        0.9 |             0.0104 |   0.0063 |    0.0148 | model-based hypothesis |
| myeloid    | C1QB         | T:cytotoxic         |        0.9 |             0.0078 |   0.0049 |    0.0112 | model-based hypothesis |
| myeloid    | C1QC         | T:cytotoxic         |        0.9 |             0.0103 |   0.0065 |    0.0144 | model-based hypothesis |
| myeloid    | CD163        | T:cytotoxic         |        0.9 |             0.0076 |   0.0044 |    0.0113 | model-based hypothesis |
| myeloid    | MRC1         | T:cytotoxic         |        0.9 |             0.0043 |   0.0022 |    0.0066 | model-based hypothesis |
| myeloid    | SELENOP      | T:cytotoxic         |        0.9 |             0.0034 |   0.0019 |    0.0051 | model-based hypothesis |
| myeloid    | SPP1         | T:cytotoxic         |        0.9 |             0.0026 |  -0.0003 |    0.0051 | model-based hypothesis |
| myeloid    | GPNMB        | T:cytotoxic         |        0.9 |             0.0015 |   0.0007 |    0.0025 | model-based hypothesis |
| myeloid    | CD9          | T:cytotoxic         |        0.9 |            -0.0014 |  -0.0022 |   -0.0006 | model-based hypothesis |
| myeloid    | FABP5        | T:cytotoxic         |        0.9 |            -0.0007 |  -0.0021 |    0.0005 | model-based hypothesis |
| myeloid    | LGALS3       | T:cytotoxic         |        0.9 |             0.0001 |  -0.0009 |    0.0009 | model-based hypothesis |

## Phase 6 pretraining probe

| receiver   | probe              | target            | R2 [CI]               | probe - strongest baseline [CI]   | probe - graph [CI]      | pretraining helps   |
|:-----------|:-------------------|:------------------|:----------------------|:----------------------------------|:------------------------|:--------------------|
| myeloid    | target_masked      | myeloid:C1Q_TAM   | 0.485 [0.382, 0.544]  | 0.151 [0.108, 0.196]              | 0.148 [0.101, 0.195]    | False               |
| myeloid    | target_masked      | myeloid:SPP1_TAM  | 0.389 [0.280, 0.466]  | 0.002 [-0.053, 0.042]             | 0.012 [-0.044, 0.046]   | False               |
| myeloid    | target_masked      | C1QA              | 0.258 [0.189, 0.299]  | 0.062 [0.051, 0.081]              | 0.066 [0.055, 0.080]    | False               |
| myeloid    | target_masked      | C1QB              | 0.334 [0.256, 0.383]  | 0.073 [0.049, 0.094]              | 0.055 [0.034, 0.076]    | False               |
| myeloid    | target_masked      | C1QC              | 0.339 [0.253, 0.400]  | 0.088 [0.074, 0.110]              | 0.075 [0.055, 0.100]    | False               |
| myeloid    | target_masked      | CD163             | 0.264 [0.166, 0.323]  | 0.075 [0.035, 0.100]              | 0.069 [0.033, 0.097]    | False               |
| myeloid    | target_masked      | MRC1              | 0.114 [0.076, 0.146]  | 0.018 [0.001, 0.032]              | 0.021 [-0.004, 0.043]   | False               |
| myeloid    | target_masked      | SELENOP           | 0.257 [0.188, 0.309]  | 0.012 [-0.024, 0.051]             | -0.014 [-0.043, 0.012]  | False               |
| myeloid    | target_masked      | SPP1              | 0.522 [0.340, 0.637]  | -0.111 [-0.173, -0.077]           | -0.132 [-0.191, -0.103] | False               |
| myeloid    | target_masked      | GPNMB             | 0.259 [0.218, 0.296]  | 0.047 [0.022, 0.075]              | 0.031 [0.006, 0.058]    | False               |
| myeloid    | target_masked      | CD9               | 0.126 [0.077, 0.156]  | 0.014 [0.001, 0.021]              | 0.001 [-0.015, 0.011]   | False               |
| myeloid    | target_masked      | FABP5             | 0.145 [0.101, 0.178]  | 0.021 [0.015, 0.030]              | 0.004 [-0.001, 0.008]   | False               |
| myeloid    | target_masked      | LGALS3            | 0.154 [0.126, 0.183]  | 0.012 [0.004, 0.025]              | 0.002 [-0.007, 0.015]   | False               |
| myeloid    | target_masked      | genes_pooled      | 0.295 [0.223, 0.345]  | 0.018 [0.003, 0.037]              | 0.005 [-0.009, 0.022]   | False               |
| myeloid    | neighbourhood_only | myeloid:C1Q_TAM   | 0.192 [0.115, 0.237]  | -0.142 [-0.196, -0.090]           | -0.144 [-0.207, -0.080] | False               |
| myeloid    | neighbourhood_only | myeloid:SPP1_TAM  | 0.240 [0.125, 0.334]  | -0.147 [-0.212, -0.107]           | -0.137 [-0.205, -0.089] | False               |
| myeloid    | neighbourhood_only | C1QA              | 0.118 [0.075, 0.156]  | -0.078 [-0.119, -0.031]           | -0.074 [-0.119, -0.032] | False               |
| myeloid    | neighbourhood_only | C1QB              | 0.182 [0.126, 0.221]  | -0.079 [-0.123, -0.035]           | -0.098 [-0.161, -0.046] | False               |
| myeloid    | neighbourhood_only | C1QC              | 0.163 [0.102, 0.200]  | -0.088 [-0.139, -0.036]           | -0.102 [-0.165, -0.045] | False               |
| myeloid    | neighbourhood_only | CD163             | 0.094 [0.027, 0.150]  | -0.096 [-0.152, -0.045]           | -0.102 [-0.165, -0.043] | False               |
| myeloid    | neighbourhood_only | MRC1              | 0.066 [0.039, 0.089]  | -0.029 [-0.039, -0.015]           | -0.026 [-0.038, -0.015] | False               |
| myeloid    | neighbourhood_only | SELENOP           | 0.119 [0.074, 0.156]  | -0.126 [-0.148, -0.088]           | -0.152 [-0.181, -0.118] | False               |
| myeloid    | neighbourhood_only | SPP1              | 0.450 [0.202, 0.581]  | -0.183 [-0.306, -0.128]           | -0.204 [-0.327, -0.145] | False               |
| myeloid    | neighbourhood_only | GPNMB             | 0.108 [0.086, 0.134]  | -0.104 [-0.147, -0.057]           | -0.120 [-0.164, -0.073] | False               |
| myeloid    | neighbourhood_only | CD9               | 0.093 [0.051, 0.120]  | -0.019 [-0.030, -0.010]           | -0.031 [-0.045, -0.022] | False               |
| myeloid    | neighbourhood_only | FABP5             | 0.109 [0.075, 0.139]  | -0.015 [-0.026, -0.002]           | -0.032 [-0.044, -0.018] | False               |
| myeloid    | neighbourhood_only | LGALS3            | 0.117 [0.094, 0.144]  | -0.025 [-0.039, -0.009]           | -0.036 [-0.049, -0.021] | False               |
| myeloid    | neighbourhood_only | genes_pooled      | 0.184 [0.120, 0.227]  | -0.093 [-0.109, -0.072]           | -0.106 [-0.126, -0.081] | False               |
| tumour     | target_masked      | tumour:classical  | 0.502 [0.342, 0.633]  | -0.192 [-0.283, -0.111]           | -0.096 [-0.172, -0.024] | False               |
| tumour     | target_masked      | tumour:basal_like | 0.359 [0.245, 0.494]  | -0.084 [-0.129, -0.044]           | 0.078 [-0.047, 0.269]   | False               |
| tumour     | target_masked      | CEACAM6           | 0.480 [0.280, 0.632]  | -0.174 [-0.370, -0.049]           | -0.172 [-0.338, -0.061] | False               |
| tumour     | target_masked      | AGR2              | 0.353 [0.247, 0.447]  | -0.108 [-0.164, -0.056]           | -0.112 [-0.151, -0.067] | False               |
| tumour     | target_masked      | S100P             | 0.379 [0.290, 0.461]  | -0.107 [-0.163, -0.059]           | -0.095 [-0.145, -0.052] | False               |
| tumour     | target_masked      | LCN2              | 0.309 [0.216, 0.396]  | -0.247 [-0.336, -0.164]           | -0.245 [-0.332, -0.164] | False               |
| tumour     | target_masked      | PIGR              | 0.389 [0.212, 0.513]  | -0.160 [-0.287, -0.076]           | -0.165 [-0.291, -0.081] | False               |
| tumour     | target_masked      | OLFM4             | 0.350 [0.013, 0.522]  | -0.292 [-0.518, -0.161]           | -0.298 [-0.530, -0.166] | False               |
| tumour     | target_masked      | SPINK1            | 0.433 [0.306, 0.541]  | -0.162 [-0.213, -0.109]           | -0.169 [-0.219, -0.118] | False               |
| tumour     | target_masked      | CEACAM1           | 0.106 [0.052, 0.155]  | -0.032 [-0.059, -0.006]           | -0.033 [-0.057, -0.003] | False               |
| tumour     | target_masked      | KRT6A             | 0.417 [0.157, 0.579]  | -0.172 [-0.219, -0.080]           | -0.173 [-0.223, -0.083] | False               |
| tumour     | target_masked      | KRT14             | 0.221 [0.103, 0.343]  | -0.086 [-0.181, -0.007]           | -0.086 [-0.187, -0.007] | False               |
| tumour     | target_masked      | KRT17             | 0.345 [0.267, 0.423]  | -0.088 [-0.147, -0.040]           | -0.068 [-0.107, -0.039] | False               |
| tumour     | target_masked      | S100A2            | 0.246 [0.050, 0.435]  | -0.210 [-0.288, -0.089]           | -0.196 [-0.249, -0.100] | False               |
| tumour     | target_masked      | SERPINB5          | 0.127 [0.056, 0.234]  | -0.026 [-0.043, -0.009]           | -0.021 [-0.033, -0.004] | False               |
| tumour     | target_masked      | AREG              | 0.058 [0.016, 0.093]  | -0.099 [-0.185, -0.039]           | -0.088 [-0.158, -0.037] | False               |
| tumour     | target_masked      | SLC2A1            | 0.218 [0.165, 0.276]  | -0.006 [-0.014, 0.005]            | -0.018 [-0.029, -0.007] | False               |
| tumour     | target_masked      | KRT15             | 0.107 [0.087, 0.128]  | 0.018 [0.008, 0.025]              | 0.024 [0.009, 0.039]    | False               |
| tumour     | target_masked      | ITGA3             | 0.263 [0.221, 0.323]  | 0.024 [0.006, 0.041]              | 0.036 [0.021, 0.053]    | False               |
| tumour     | target_masked      | MET               | 0.074 [0.046, 0.100]  | -0.016 [-0.034, -0.003]           | -0.002 [-0.012, 0.007]  | False               |
| tumour     | target_masked      | genes_pooled      | 0.310 [0.246, 0.386]  | -0.125 [-0.158, -0.092]           | -0.125 [-0.156, -0.096] | False               |
| tumour     | neighbourhood_only | tumour:classical  | 0.311 [0.217, 0.390]  | -0.383 [-0.474, -0.283]           | -0.286 [-0.401, -0.150] | False               |
| tumour     | neighbourhood_only | tumour:basal_like | 0.177 [0.036, 0.323]  | -0.266 [-0.354, -0.195]           | -0.104 [-0.277, 0.109]  | False               |
| tumour     | neighbourhood_only | CEACAM6           | 0.284 [0.180, 0.384]  | -0.369 [-0.475, -0.237]           | -0.367 [-0.476, -0.242] | False               |
| tumour     | neighbourhood_only | AGR2              | 0.256 [0.183, 0.333]  | -0.205 [-0.269, -0.126]           | -0.209 [-0.269, -0.132] | False               |
| tumour     | neighbourhood_only | S100P             | 0.272 [0.177, 0.371]  | -0.214 [-0.274, -0.158]           | -0.203 [-0.264, -0.152] | False               |
| tumour     | neighbourhood_only | LCN2              | 0.209 [0.122, 0.303]  | -0.347 [-0.435, -0.245]           | -0.345 [-0.432, -0.244] | False               |
| tumour     | neighbourhood_only | PIGR              | 0.188 [0.093, 0.274]  | -0.361 [-0.460, -0.249]           | -0.366 [-0.466, -0.252] | False               |
| tumour     | neighbourhood_only | OLFM4             | 0.106 [-0.052, 0.267] | -0.536 [-0.706, -0.339]           | -0.541 [-0.706, -0.355] | False               |
| tumour     | neighbourhood_only | SPINK1            | 0.297 [0.160, 0.444]  | -0.298 [-0.430, -0.166]           | -0.305 [-0.439, -0.175] | False               |
| tumour     | neighbourhood_only | CEACAM1           | 0.039 [-0.046, 0.097] | -0.099 [-0.160, -0.053]           | -0.100 [-0.154, -0.053] | False               |
| tumour     | neighbourhood_only | KRT6A             | 0.213 [0.074, 0.297]  | -0.375 [-0.501, -0.182]           | -0.377 [-0.507, -0.170] | False               |
| tumour     | neighbourhood_only | KRT14             | 0.144 [0.056, 0.247]  | -0.163 [-0.278, -0.041]           | -0.162 [-0.284, -0.040] | False               |
| tumour     | neighbourhood_only | KRT17             | 0.240 [0.144, 0.363]  | -0.194 [-0.278, -0.112]           | -0.173 [-0.242, -0.107] | False               |
| tumour     | neighbourhood_only | S100A2            | 0.111 [-0.086, 0.257] | -0.345 [-0.433, -0.179]           | -0.331 [-0.442, -0.177] | False               |
| tumour     | neighbourhood_only | SERPINB5          | 0.105 [0.064, 0.159]  | -0.048 [-0.108, -0.004]           | -0.043 [-0.102, 0.001]  | False               |
| tumour     | neighbourhood_only | AREG              | 0.057 [-0.017, 0.120] | -0.101 [-0.151, -0.046]           | -0.090 [-0.131, -0.043] | False               |
| tumour     | neighbourhood_only | SLC2A1            | 0.157 [0.117, 0.204]  | -0.067 [-0.088, -0.042]           | -0.079 [-0.105, -0.049] | False               |
| tumour     | neighbourhood_only | KRT15             | 0.078 [0.066, 0.091]  | -0.011 [-0.027, 0.001]            | -0.005 [-0.026, 0.015]  | False               |
| tumour     | neighbourhood_only | ITGA3             | 0.183 [0.131, 0.252]  | -0.057 [-0.092, -0.022]           | -0.044 [-0.077, -0.008] | False               |
| tumour     | neighbourhood_only | MET               | 0.067 [0.043, 0.091]  | -0.023 [-0.037, -0.008]           | -0.009 [-0.017, -0.001] | False               |
| tumour     | neighbourhood_only | genes_pooled      | 0.188 [0.150, 0.234]  | -0.247 [-0.305, -0.180]           | -0.247 [-0.305, -0.183] | False               |

## Phase 6 learning curve

| receiver   | target           |   n_train_tumours |   mean_train_cells |   mean_r2 |   n_folds |   loglinear_slope_per_decade |   r2_at_10x_cells_extrapolated | note                                        |
|:-----------|:-----------------|------------------:|-------------------:|----------:|----------:|-----------------------------:|-------------------------------:|:--------------------------------------------|
| myeloid    | myeloid:C1Q_TAM  |                 3 |            21517   |    0.2615 |         5 |                       0.0639 |                         0.3628 | descriptive log-linear extrapolation; rough |
| myeloid    | myeloid:C1Q_TAM  |                 6 |            47024.2 |    0.2983 |         5 |                       0.0639 |                         0.3628 | descriptive log-linear extrapolation; rough |
| myeloid    | myeloid:C1Q_TAM  |                 7 |            54330   |    0.2002 |         1 |                       0.0639 |                         0.3628 | descriptive log-linear extrapolation; rough |
| myeloid    | myeloid:C1Q_TAM  |                 8 |            64082   |    0.2853 |         1 |                       0.0639 |                         0.3628 | descriptive log-linear extrapolation; rough |
| myeloid    | myeloid:C1Q_TAM  |                 9 |            52765   |    0.4002 |         2 |                       0.0639 |                         0.3628 | descriptive log-linear extrapolation; rough |
| tumour     | tumour:classical |                 3 |            22802.4 |    0.4931 |         5 |                       0.066  |                         0.5879 | descriptive log-linear extrapolation; rough |
| tumour     | tumour:classical |                 6 |            47730.4 |    0.52   |         5 |                       0.066  |                         0.5879 | descriptive log-linear extrapolation; rough |
| tumour     | tumour:classical |                 7 |            60778   |    0.4508 |         1 |                       0.066  |                         0.5879 | descriptive log-linear extrapolation; rough |
| tumour     | tumour:classical |                 8 |            57906   |    0.5298 |         1 |                       0.066  |                         0.5879 | descriptive log-linear extrapolation; rough |
| tumour     | tumour:classical |                 9 |            64704.5 |    0.5725 |         2 |                       0.066  |                         0.5879 | descriptive log-linear extrapolation; rough |
