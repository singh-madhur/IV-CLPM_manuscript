# IV-CLPM_manuscript

This is a companion repository for the code used for simulations and analyses in the manuscript:


Singh, M., Verhulst, B., Vinh, P., Zhou, Y., Castro-de-Araujo, L. F. S., Maes, H. H., Dolan, C. V., and Neale, M. C. (Under Review). Integrating Cross-Lagged Panel Models with Instrumental Variables to Extend the Temporal Generalizability of Causal Inference.


Abstract: 

Cross-lagged panel models (CLPMs) are commonly used to estimate causal influences between two variables with repeated measurements. In a CLPM, the magnitude and significance of the lagged effects depend on the time interval between measurement occasions, and these effects usually become undetectable at longer intervals. To address this limitation, we integrated instrumental variables and the CLPM, allowing us to estimate simultaneously the lagged (i.e., “distal”) effects and the bidirectional cross-sectional (i.e., “proximal”) effects at each wave. The distal effects reflect Granger causal effects across time, which decay with increasing time intervals. The proximal effects at each wave capture causal influences that have accrued over time and can help infer causality when the distal effects become undetectable at longer intervals. The distal and proximal causal estimates can be tested separately or jointly through likelihood-ratio tests, unraveling distinct temporal aspects of the causal process. Significant proximal effects, with negligible distal effects, imply that the time interval is too long for studying Granger causality in a pair of variables. We demonstrate that the time interval between measurement occasions profoundly impacts the ability to identify causal effects in longitudinal data and propose strategies to detect causation regardless of the time interval in a study. 


The directory `simulations` contains the scripts for simulating data and fitting the models. The directory `manuscript_figures_n_tables` contains the scripts for the figures and table presented in the manuscript.


For queries, please contact Madhur Singh at "singhm18 at vcu dot edu".
