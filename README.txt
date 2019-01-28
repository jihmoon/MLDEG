MLDEG v1.0	by Ji Hwan Moon

1. netFeatures.py
	python netFeatures.py -n <STRING PPI in entrezID> -i <Gene Expression Profile in entrezid> -r <PCC cutoff> -p <PCC p-value cutoff> -o <Output File Prefix> -s <# of treated samples:# of control samples>
	Ex) python netFeatures.py -n 9606.protein.entrez.v10.5.txt -i GSE90116.profile -r 0.7 -p 0.05 -o GSE90116 -s 3:3
	- Input Files:
		- STRING PPI in entrezID:
			"
			SOURCE	TARGET	SCORE
			...
			"
			- Ex)
				9997    94241   281
				9997    9601    371
				9997    9617    183
				9997    9636    247
				9997    9669    153
				9997    9804    154
				9997    9843    478
				9997    9940    189
				9997    995     229
				9997    9973    187
				...

		- Gene Expression Profile Format:
			"
				Treated_replicate1	...	Treated_replicateM	Control_replicate1	...	Control_replicateN
			EntrezID	EXP_VALUE	...	EXP_VALUE	EXP_VALUE	...	EXP_VALUE
			...
			"
			- Ex)
				        10th_s9wt       6th_s11wt       16th_s14wt      10th_s7ko       15th_s12ko      16th_s12ko
				11287   1.379933488     2.747453629     0.540946188     1.228285999     0.235888783     1.363161581
				11298   0       0       0       0       0       0
				11302   7.341523198     10.53064513     7.389169431     7.804707203     8.656584081     10.91453295
				11303   25.93261111     21.9361143      27.39481197     12.88898202     7.210632674     10.40312353
				11304   0       0       0       0.13221134      0.152344839     0
				11305   3.853984408     2.23361402      5.261609296     8.291673316     4.246378506     4.461655787
				11306   5.23427454      2.456487842     3.407589959     6.156608502     8.244556666     3.726901708
				11307   5.657018011     0.442104862     6.854875117     4.98239958      6.832422499     5.982336347
				11308   12.64083886     29.7296594      19.15025152     9.093151865     22.37382215     22.34925056
				...

		- Output Files:
		- [Output Prefix]_corr.stats.txt
		- [Output Prefix]_node.corr.txt
		- [Output Prefix]_weighted.ppi
		- [Output Prefix]_netProp.txt
		- Ex)
			- GSE90116_corr.stats.txt
			- GSE90116_node.corr.txt
			- GSE90116_weighted.ppi
			- GSE90116_netProp.txt

2. MLDEG.pl
	perl MLDEG.pl <DEG Results> <PCC Result> <Network Propagation Result> <log2FC cutoff for training data> <p-value cutoff for training data> <Output File Prefix>
	Ex) perl MLDEG.pl GSE90116_deg.results GSE90116_corr.stats.txt GSE90116_netProp.txt 0.5 0.01 GSE90116

	- Input Files:
		- DEG Results Format:
			"
			EntrezID	p-value_Tool1	log2FC_Tool1	p-value_Tool2	log2FC_Tool2	...	p-value_ToolK	log2FC_ToolK
			...
			"
			- Ex)
			11287   0.941076821654928       0.477039638273295       1       1.36460536947671        0.999999999986711       0.564133080400117       0.736750355958965       0.65947320088546
			11302   0.952709324281507       -0.106845234063094      1       -0.122685393751474      0.999999999995944       -0.0971657518516668     0.908874208666597       -0.108426866533915
			11303   0.272982279609671       1.3019751884409 1       1.36018665792367        0.999999945984339       1.28506391572805        0.0103738408139929      1.29899103570804
			11304   0.689373695998959       -0.235616743996669      1       0       0.999999999886136       -0.325724666362309      1       -0.813133293616092
			...
		- PCC Result:
			- [Output Prefix]_corr.stats.txt of netFeatures.py
		- Network Propagation Result
			- [Output Prefix]_netProp.txt of netFeatures.py

	- Output Files:
		- [Output Prefix]_trainingset.list
		- [Output Prefix]_testset.list
		- [Output Prefix]_TR.arff
		- [Output Prefix]_TS.arff
		- Ex)
			- GSE90116_trainingset.list
			- GSE90116_testset.list
			- GSE90116_TR.arff
			- GSE90116_TS.arff

3. weka.runner.pl
	perl weka.runner.pl <Classifier (SMO, LR or J48)> <Output File Prefix>
	Ex) perl weka.runner.pl LR GSE90116
	- Output Files:
		- [Output Prefix].[Classifier].TR.results
		- [Output Prefix].[Classifier].results
		- [Output Prefix].[Classifier].model
		- Ex)
			- GSE90116.LR.TR.results
			- GSE90116.LR.results
			- GSE90116.LR.model

4. weka.prediction.parser.pl
	perl weka.prediction.parser.pl <Test Gene Set> <Weka Prediction Result> <Output File Prefix> <Confidence Cutoff>
	Ex) perl weka.prediction.parser.pl GSE90116_testset.list GSE90116.LR.results GSE90116 0
	- Output Files:
		- [Output Prefix]_predicted.result
		- Ex)
			- GSE90116_predicted.result
5. MLDEG.results.pl
	perl MLDEG.results.pl <Training Gene Set> <Predicted Gene Set> <Output File Prefix>
	Ex) perl MLDEG.results.pl GSE90116_trainingset.list GSE90116_predicted.result GSE90116
	- Output Files:
		- [Output Prefix].final.result
		- Ex)
			GSE90116.final.result

