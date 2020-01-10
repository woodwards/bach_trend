# prepare output dataframes

#
parcols <- c('meda0raw', 'medb0raw', 'medc1raw', 'medc2raw', 'medd1raw', 'medd2raw', 
             'slowa0raw', 'slowb0raw', 'medc1raw', 'medc2raw', 'slowd1raw', 'medd2raw', 
             'chem1i', 'chem1wt', 'chem1ae', 'chem1re', 
             'chem1fastraw0', 'chem1medraw0', 'chem1slowraw0', 'chem1fastraw1', 'chem1medraw1', 'chem1slowraw1',
             'chem1fastrawb1', 'chem1medrawb1', 'chem1slowrawb1', 'chem1fastrawb2', 'chem1medrawb2', 'chem1slowrawb2',
             'chem2i', 'chem2wt', 'chem2ae', 'chem2re', 
             'chem2fastraw0', 'chem2medraw0', 'chem2slowraw0', 'chem2fastraw1', 'chem2medraw1', 'chem2slowraw1',
             'chem2fastrawb1', 'chem2medrawb1', 'chem2slowrawb1', 'chem2fastrawb2', 'chem2medrawb2', 'chem2slowrawb2',
             'medth', 'slowth', 'spare', 'spare', 'spare', 'area')
statcols <- c('llchem1', 'llchem2', 'lltotal', 'vllchem1', 'vllchem2', 'vlltotal',
              'fastflow', 'medflow', 'slowflow', 'totalflow', 'fastflowpc', 'medflowpc', 'slowflowpc', 
              'fastdays', 'meddays', 'slowdays',
              'nse1', 'nse2', 'vnse1', 'vnse2',
              'rmse1', 'rmse2', 'vrmse1', 'vrmse2',
              'gme1', 'gme2', 'vgme1', 'vgme2',
              'grmse1', 'grmse2', 'vgrmse1', 'vgrmse2',
              'meda0', 'meda1', 'meda2', 'medb0', 'medb1', 'medb2',
              'slowa0', 'slowa1', 'slowa2', 'slowb0', 'slowb1', 'slowb2',
              'medBFImax', 'slowBFImax', 'medrec', 'slowrec',
              'chem1fast0', 'chem1med0', 'chem1slow0',
              'chem2fast0', 'chem2med0', 'chem2slow0',
              'chem1fast1', 'chem1med1', 'chem1slow1',
              'chem2fast1', 'chem2med1', 'chem2slow1',
              'chem1fastb1', 'chem1medb1', 'chem1slowb1',
              'chem2fastb1', 'chem2medb1', 'chem2slowb1',
              'chem1fastb2', 'chem1medb2', 'chem1slowb2',
              'chem2fastb2', 'chem2medb2', 'chem2slowb2',
              'totalflow', 
              'totalgens', 'elapsed', 'gelmanr',
              'nTP', 'nTN', 'vnTP', 'vnTN',
              'fastTPload', 'medTPload', 'slowTPload', 'totalTPload',
              'fastTNload', 'medTNload', 'slowTNload', 'totalTNload',
              'fastTPw', 'medTPw', 'slowTPw', 'totalTPw',
              'fastTNw', 'medTNw', 'slowTNw', 'totalTNw',
              'fastTPloadpc', 'medTPloadpc', 'slowTPloadpc', 
              'fastTNloadpc', 'medTNloadpc', 'slowTNloadpc', 
              'auto1', 'auto2', 'vauto1', 'vauto2')
if (from_scratch) {
  bestxs <- as.data.frame(matrix(data=0, ncol = length(parcols), nrow = nruns, dimnames=list(NULL,parcols)))
  quartxs <- as.data.frame(matrix(data=0, ncol = length(parcols), nrow = nruns, dimnames=list(NULL,parcols)))
  sbestxs <- as.data.frame(matrix(data=0, ncol = length(statcols), nrow = nruns, dimnames=list(NULL,statcols)))
  squartxs <- as.data.frame(matrix(data=0, ncol = length(statcols), nrow = nruns, dimnames=list(NULL,statcols)))
  runrecord <- runlist
  runrecord$gelmanr <- NA
  runrecord$elapsed <- NA
} else {  
  bestxs <- read_tsv(paste(out_path, 'bestxs.tsv', sep=''), col_names=FALSE, col_types=cols())
  quartxs <- read_tsv(paste(out_path, 'quartxs.tsv', sep=''), col_names=FALSE, col_types=cols())
  sbestxs <- read_tsv(paste(out_path, 'sbestxs.tsv', sep=''), col_names=FALSE, col_types=cols())
  squartxs <- read_tsv(paste(out_path, 'squartxs.tsv', sep=''), col_names=FALSE, col_types=cols())
  runrecord <- read_tsv(paste(out_path, 'runrecord.tsv', sep=''), col_names=TRUE, col_types=cols())
}

