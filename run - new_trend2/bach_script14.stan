functions { // internally defined
  
}

data { // externally supplied

  int nintoptions;                 // n integer options
  int intoptions[nintoptions];     // integer options (dates)
  int nrealoptions;                // n real options
  real realoptions[nrealoptions];  // real options (errors)
  int ndate;                       // number of dates  
  int date[ndate];                 // dates
  real flow[ndate];                 // flows
  real TP[ndate];                   // TP data (NA is coded as -1)
  real TN[ndate];                   // TN data (NA is coded as -1)
  real meanTPcalib;
  real meanTNcalib;

} 

transformed data { // internally generated

  // unpack intoptions here since integers can't be considered parameters
  int startrun = intoptions[1];
  int startcalib = intoptions[2];
  int endcalib = intoptions[3];
  int startvalid = intoptions[4];
  int endvalid = intoptions[5];
  
  // useful value
  real llzero = student_t_lpdf(0 | 7, 0, 1);
  real eps = machine_precision();
  
}

parameters { // raw parameters to be adjusted

  // declare raw adjustable parameters here
  // good idea if these all on similar scale
  
  real<lower=0.0, upper=1.0> medb0raw;
  real<lower=0.0, upper=1.0> medd1raw;       
  real<lower=0.0, upper=1.0> slowb0raw;
  real<lower=0.0, upper=1.0> slowd1raw;       
  real<lower=0.0, upper=1.0> chem1fastraw0; // initial conc     
  real<lower=0.0, upper=1.0> chem1medraw0;
  real<lower=0.0, upper=1.0> chem1slowraw0;
  real<lower=0.0, upper=1.0> chem2fastraw0;
  real<lower=0.0, upper=1.0> chem2medraw0;
  real<lower=0.0, upper=1.0> chem2slowraw0;
  real<lower=0.0, upper=1.0> chem1fastraw1; // final conc  
  real<lower=0.0, upper=1.0> chem1medraw1;
  real<lower=0.0, upper=1.0> chem1slowraw1;
  real<lower=0.0, upper=1.0> chem2fastraw1;
  real<lower=0.0, upper=1.0> chem2medraw1;
  real<lower=0.0, upper=1.0> chem2slowraw1;
  
  real<lower=-1.0, upper=1.0> chem1fastrawb1; // first harmonic
  real<lower=-1.0, upper=1.0> chem1medrawb1;
  real<lower=-1.0, upper=1.0> chem1slowrawb1;
  real<lower=-1.0, upper=1.0> chem2fastrawb1;
  real<lower=-1.0, upper=1.0> chem2medrawb1;
  real<lower=-1.0, upper=1.0> chem2slowrawb1;
  real<lower=-1.0, upper=1.0> chem1fastrawb2; // second harmonic 
  real<lower=-1.0, upper=1.0> chem1medrawb2;
  real<lower=-1.0, upper=1.0> chem1slowrawb2;
  real<lower=-1.0, upper=1.0> chem2fastrawb2;
  real<lower=-1.0, upper=1.0> chem2medrawb2;
  real<lower=-1.0, upper=1.0> chem2slowrawb2;
  
} 

transformed parameters { // additional parameters to be used by the model

  // put fixed parameters here
  // real medb0raw = 1;
  // real medd1raw = 0;       
  // real slowb0raw = 1;
  // real chem1fastraw = 0;     
  // real chem2fastraw = 0;
  
  // unpack these here so available in output
  real chem1ae = realoptions[1];
  real chem1re = realoptions[2];
  real chem2ae = realoptions[3];
  real chem2re = realoptions[4];
  real area = realoptions[5];
  
  real priorwide = 0.3;
  real priornarrow = 0.1;
  
  real medd1scale;
  real meda1;
  real medb0;
  real medBFImax;
  real medrec;
  real slowd1scale;
  real slowa1;
  real slowb0;
  real slowBFImax;
  real slowrec;
  real ratio;
  real bach[ndate,11]; // need this here to calc auxiliaries later
  real alpha;

  real chem1fast0;     
  real chem1med0;
  real chem1slow0;
  real chem2fast0;
  real chem2med0;
  real chem2slow0;
  real chem1fast1;     
  real chem1med1;
  real chem1slow1;
  real chem2fast1;
  real chem2med1;
  real chem2slow1;

  real chem1fastb1;     
  real chem1medb1;
  real chem1slowb1;
  real chem2fastb1;
  real chem2medb1;
  real chem2slowb1;
  real chem1fastb2;     
  real chem1medb2;
  real chem1slowb2;
  real chem2fastb2;
  real chem2medb2;
  real chem2slowb2;
  
  // real<lower=0.5, upper=1.0> medb0raw;
  // real<lower=-0.1*log(1-0.5), upper=-0.1*log(1-0.99)> medd1raw;       
  // real<lower=0.0, upper=1.0> slowb0raw;
  // real<lower=-0.1*log(1-0.99), upper=-0.1*log(1-0.9999)> slowd1raw;       
  
  medd1scale = (-0.1*log(1-0.5)) + ((-0.1*log(1-0.99)) - (-0.1*log(1-0.5))) * medd1raw; // = 0.1*ln(1-meda1)
  meda1 = (1-exp(-fabs(medd1scale)*10)); // for d1>0
  ratio = 1/(1-meda1);
  medb0 = ratio>1 ? medb0raw/ratio : medb0raw;
  medBFImax = medb0/(1-meda1);
  medrec = medb0<1 ? meda1/(1-medb0) : 0;
  
  slowd1scale = (-0.1*log(1-0.99)) + ((-0.1*log(1-0.9999)) - (-0.1*log(1-0.99))) * slowd1raw;// = 0.1*ln(1-slowa1)
  slowa1 = (1-exp(-fabs(slowd1scale)*10)); // for d1>0
  ratio = 1/(1-slowa1);
  slowb0 = ratio>1 ? slowb0raw/ratio : slowb0raw;
  slowBFImax = slowb0/(1-slowa1);
  slowrec = slowb0<1 ? slowa1/(1-slowb0) : 0;

  chem1fast0 = chem1fastraw0 * 2;
  chem1med0 = chem1medraw0 * 2;
  chem1slow0 = chem1slowraw0 * 2;
  chem2fast0 = chem2fastraw0 * 12;
  chem2med0 = chem2medraw0 * 12;
  chem2slow0 = chem2slowraw0 * 12;
  
  chem1fast1 = chem1fastraw1 * 2;
  chem1med1 = chem1medraw1 * 2;
  chem1slow1 = chem1slowraw1 * 2;
  chem2fast1 = chem2fastraw1 * 12;
  chem2med1 = chem2medraw1 * 12;
  chem2slow1 = chem2slowraw1 * 12;

  chem1fastb1 = chem1fastrawb1 * 2;
  chem1medb1 = chem1medrawb1 * 2;
  chem1slowb1 = chem1slowrawb1 * 2;
  chem2fastb1 = chem2fastrawb1 * 12;
  chem2medb1 = chem2medrawb1 * 12;
  chem2slowb1 = chem2slowrawb1 * 12;

  chem1fastb2 = chem1fastrawb2 * 2;
  chem1medb2 = chem1medrawb2 * 2;
  chem1slowb2 = chem1slowrawb2 * 2;
  chem2fastb2 = chem2fastrawb2 * 12;
  chem2medb2 = chem2medrawb2 * 12;
  chem2slowb2 = chem2slowrawb2 * 12;
  
  // initialise model output array
  for (n in 1:ndate) 
    for (i in 1:11) 
      bach[n,i] = 0;
     
  // initial conditions
  bach[startrun,1] = flow[startrun];
  bach[startrun,2] = flow[startrun]*medBFImax;
  bach[startrun,3] = flow[startrun]*medBFImax*slowBFImax;
  bach[startrun+1,1] = flow[startrun+1];
  bach[startrun+1,2] = flow[startrun+1]*medBFImax;
  bach[startrun+1,3] = flow[startrun+1]*medBFImax*slowBFImax;

  for (n in (startrun+2):endvalid) {
    
    // double flow filter
    bach[n,1] = flow[n];
    bach[n,2] = meda1*bach[n-1,2] + medb0*bach[n,1];
    if (bach[n,2] > bach[n,1]) 
      bach[n,2] = bach[n,1];
    bach[n,3] = slowa1*bach[n-1,3] + slowb0*bach[n,2];
    if (bach[n,3] > bach[n,2]) 
      bach[n,3] = bach[n,2];
    
    // flow path concs
    alpha = (n-366.0)/(5844.0-366.0);
    bach[n,6] = fmax(eps, chem1fast0*(1-alpha)+chem1fast1*alpha+chem1fastb1*sin(pi()*alpha)+chem1fastb2*sin(2*pi()*alpha));
    bach[n,7] = fmax(eps, chem1med0 *(1-alpha)+chem1med1 *alpha+chem1medb1*sin(pi()*alpha)+chem1medb2*sin(2*pi()*alpha));
    bach[n,8] = fmax(eps, chem1slow0*(1-alpha)+chem1slow1*alpha+chem1slowb1*sin(pi()*alpha)+chem1slowb2*sin(2*pi()*alpha));
    bach[n,9] = fmax(eps, chem2fast0*(1-alpha)+chem2fast1*alpha+chem2fastb1*sin(pi()*alpha)+chem2fastb2*sin(2*pi()*alpha));
    bach[n,10] = fmax(eps, chem2med0 *(1-alpha)+chem2med1 *alpha+chem2medb1*sin(pi()*alpha)+chem2medb2*sin(2*pi()*alpha));
    bach[n,11] = fmax(eps, chem2slow0*(1-alpha)+chem2slow1*alpha+chem2slowb1*sin(pi()*alpha)+chem2slowb2*sin(2*pi()*alpha));
    
    // overall concs
    bach[n,4] = ((bach[n,1]-bach[n,2])*bach[n,6] + 
                 (bach[n,2]-bach[n,3])*bach[n,7] + 
                  bach[n,3]           *bach[n,8] ) / bach[n,1];
    bach[n,5] = ((bach[n,1]-bach[n,2])*bach[n,9] + 
                 (bach[n,2]-bach[n,3])*bach[n,10] + 
                  bach[n,3]           *bach[n,11] ) / bach[n,1];
    
  } // model loop

}

model { // compare model with data

  // variables declared here are local
  real s;
  real z;
  
  // priors
  medb0raw ~ normal(1.0, priornarrow); // see Priors.xlsx for calculations
  medd1raw ~ normal(0.4, priorwide); 
  slowb0raw ~ normal(1.0, priornarrow);
  slowd1raw ~ normal(0.7, priorwide); 
  chem1fastraw0 ~ normal(0.2/2.0, priorwide); 
  chem1medraw0 ~ normal(0.1/2.0, priorwide); 
  chem1slowraw0 ~ normal(0.1/2.0, priorwide); 
  chem2fastraw0 ~ normal(2.0/12.0, priorwide);
  chem2medraw0 ~ normal(3.0/12.0, priorwide); 
  chem2slowraw0 ~ normal(1.0/12.0, priorwide);
  chem1fastraw1 ~ normal(0.2/2.0, priorwide); 
  chem1medraw1 ~ normal(0.1/2.0, priorwide); 
  chem1slowraw1 ~ normal(0.1/2.0, priorwide); 
  chem2fastraw1 ~ normal(2.0/12.0, priorwide);
  chem2medraw1 ~ normal(3.0/12.0, priorwide); 
  chem2slowraw1 ~ normal(1.0/12.0, priorwide);
  chem1fastrawb1 ~ normal(0, priornarrow); 
  chem1medrawb1 ~ normal(0, priornarrow); 
  chem1slowrawb1 ~ normal(0, priornarrow); 
  chem2fastrawb1 ~ normal(0, priornarrow);
  chem2medrawb1 ~ normal(0, priornarrow); 
  chem2slowrawb1 ~ normal(0, priornarrow);
  chem1fastrawb2 ~ normal(0, priornarrow); 
  chem1medrawb2 ~ normal(0, priornarrow); 
  chem1slowrawb2 ~ normal(0, priornarrow); 
  chem2fastrawb2 ~ normal(0, priornarrow);
  chem2medrawb2 ~ normal(0, priornarrow); 
  chem2slowrawb2 ~ normal(0, priornarrow);
  
  // likelihood
  for (n in startcalib:endcalib) { // fit definition

    if (TP[n] >= 0) {
      s = chem1ae + chem1re*bach[n,4];
      z = (TP[n] - bach[n,4]) / s;
      z ~ student_t(7, 0, 1); 
//      TP[n] ~ student_t(7, bach[n,4], s); 
    }
    
    if (TN[n] >= 0) {
      s = chem2ae + chem2re*bach[n,5];
      z = (TN[n] - bach[n,5]) / s;
      z ~ student_t(7, 0, 1); 
//      TN[n] ~ student_t(7, bach[n,5], s); 
    }

  } // fit definition
  
} 

generated quantities { // additional values of interest

  int nTP = 0;
  int nTN = 0;
  int vnTP = 0;
  int vnTN = 0;
  int totaldays = endcalib - startcalib + 1;
  real fastflow = 0;
  real medflow = 0;
  real slowflow = 0;
  real totalflow;
  real fastflowpc;
  real medflowpc;
  real slowflowpc;
  real fastTPload = 0;
  real medTPload = 0;
  real slowTPload = 0;
  real totalTPload;
  real fastTPw = 0;
  real medTPw = 0;
  real slowTPw = 0;
  real totalTPw;
  real fastTPloadpc;
  real medTPloadpc;
  real slowTPloadpc;
  real fastTNload = 0;
  real medTNload = 0;
  real slowTNload = 0;
  real totalTNload;
  real fastTNw = 0;
  real medTNw = 0;
  real slowTNw = 0;
  real totalTNw;
  real fastTNloadpc;
  real medTNloadpc;
  real slowTNloadpc;
  real llchem1 = 0; //
  real llchem2 = 0;
  real lltotal;
  real vllchem1 = 0; 
  real vllchem2 = 0;
  real vlltotal;
  real sllchem1 = 0; // less llzero
  real sllchem2 = 0;
  real vsllchem1 = 0; // less llzero
  real vsllchem2 = 0;
  real grmse1; 
  real grmse2;
  real vgrmse1; 
  real vgrmse2;
  real auto1 = 0; 
  real auto2 = 0;
  real vauto1 = 0; 
  real vauto2 = 0;
  int totalgens = -1; // need these to replace later
  int gelmann = -1;
  real gelmanr = -1;

  real s;
  
  real z1;
  real lastz1;
  real meanz1;
  real vmeanz1;
  real llz1;
  real sumz1 = 0;
  real sumz12 = 0;
  real sumz1a = 0;
  real sumz1b = 0;
  real vsumz1 = 0;
  real vsumz12 = 0;
  real vsumz1a = 0;
  real vsumz1b = 0;

  real z2;
  real lastz2;
  real meanz2;
  real vmeanz2;
  real llz2;
  real sumz2 = 0;
  real sumz22 = 0;
  real sumz2a = 0;
  real sumz2b = 0;
  real vsumz2 = 0;
  real vsumz22 = 0;
  real vsumz2a = 0;
  real vsumz2b = 0;
  
  // do calibration calcs
  lastz1 = not_a_number();
  lastz2 = not_a_number();
  for (n in startcalib:endcalib) {

    fastflow = fastflow + (bach[n,1] - bach[n,2]); // m3/s
    medflow = medflow + (bach[n,2] - bach[n,3]);
    slowflow = slowflow + bach[n,3]; 
    
    fastTPload = fastTPload + (bach[n,1] - bach[n,2]) * bach[n,6]; // m3/s * mg/L
    medTPload = medTPload + (bach[n,2] - bach[n,3]) * bach[n,7];
    slowTPload = slowTPload + bach[n,3] * bach[n,8];

    fastTNload = fastTNload + (bach[n,1] - bach[n,2]) * bach[n,9];
    medTNload = medTNload + (bach[n,2] - bach[n,3]) * bach[n,10];
    slowTNload = slowTNload + bach[n,3] * bach[n,11];

    if (TP[n] >= 0) {
      nTP = nTP + 1;
      s = chem1ae + chem1re*bach[n,4];
      z1 = (TP[n] - bach[n,4]) / s;
      sumz1 = sumz1 + z1;
      sumz12 = sumz12 + z1 * z1;
      if (!is_nan(lastz1)) {
        sumz1a = sumz1a + z1 * lastz1;
        sumz1b = sumz1b + z1 + lastz1;
      }
      lastz1 = z1;
      llz1 = student_t_lpdf(z1 | 7, 0, 1);
      llchem1 = llchem1 + llz1; 
      sllchem1 = sllchem1 + 2 * s ^ 2 * (llzero - llz1);
    }

    if (TN[n] >= 0) {
      nTN = nTN + 1;
      s = chem2ae + chem2re*bach[n,5];
      z2 = (TN[n] - bach[n,5]) / s;
      sumz2 = sumz2 + z2;
      sumz22 = sumz22 + z2 * z2;
      if (!is_nan(lastz2)) {
        sumz2a = sumz2a + z2 * lastz2;
        sumz2b = sumz2b + z2 + lastz2;
      }
      lastz2 = z2;
      llz2 = student_t_lpdf(z2 | 7, 0, 1);
      llchem2 = llchem2 + llz2; 
      sllchem2 = sllchem2 + 2 * s ^ 2 * (llzero - llz2); 
    }
    
  } // calibration calcs
  
  meanz1 = nTP>2 ? sumz1/nTP : 0;
  meanz2 = nTN>2 ? sumz2/nTN : 0;
  auto1 = nTP>2 ? (sumz1a-meanz1*sumz1b+(nTP-1)*meanz1^2)/(sumz12-nTP*meanz1^2) : 0;
  auto2 = nTN>2 ? (sumz2a-meanz2*sumz2b+(nTN-1)*meanz2^2)/(sumz22-nTN*meanz2^2) : 0;

  fastflow = fastflow/totaldays; // m3/s
  medflow = medflow/totaldays;
  slowflow = slowflow/totaldays;
  totalflow = fastflow + medflow + slowflow;
  fastflowpc = fastflow/totalflow;
  medflowpc = medflow/totalflow;
  slowflowpc = slowflow/totalflow;
  
  fastTPload = fastTPload * 0.0864; // t
  medTPload = medTPload * 0.0864; // t
  slowTPload = slowTPload * 0.0864; // t
  totalTPload = fastTPload + medTPload + slowTPload;
  fastTPloadpc = fastTPload / totalTPload;
  medTPloadpc = medTPload / totalTPload;
  slowTPloadpc = slowTPload / totalTPload;
  
  fastTPw = 0; // g/m3
  medTPw = 0; // g/m3
  slowTPw = 0; // g/m3
  totalTPw = fastTPw + medTPw + slowTPw; // g/m3
  
  fastTNload = fastTNload * 0.0864; // t
  medTNload = medTNload * 0.0864; // t
  slowTNload = slowTNload * 0.0864; // t
  totalTNload = fastTNload + medTNload + slowTNload;
  fastTNloadpc = fastTNload / totalTNload;
  medTNloadpc = medTNload / totalTNload;
  slowTNloadpc = slowTNload / totalTNload;

  fastTNw = 0; // g/m3
  medTNw = 0; // g/m3
  slowTNw = 0; // g/m3
  totalTNw = fastTNw + medTNw + slowTNw; // g/m3
  
  llchem1 = -llchem1;
  llchem2 = -llchem2;
  lltotal = llchem1 + llchem2;
  grmse1 = nTP>2 ? sqrt(sllchem1/nTP) : 0;
  grmse2 = nTN>2 ? sqrt(sllchem2/nTN) : 0;
  
  // do validation calcs
  lastz1 = not_a_number();
  lastz2 = not_a_number();
  for (n in startvalid:endvalid) {

    if (TP[n] >= 0) {
      vnTP  = vnTP + 1;
      s = chem1ae + chem1re*bach[n,4];
      z1 = (TP[n] - bach[n,4]) / s;
      vsumz1 = vsumz1 + z1;
      vsumz12 = vsumz12 + z1 * z1;
      if (!is_nan(lastz1)) {
        vsumz1a = vsumz1a + z1 * lastz1;
        vsumz1b = vsumz1b + z1 + lastz1;
      }
      lastz1 = z1;
      llz1 = student_t_lpdf(z1 | 7, 0, 1);
      vllchem1 = vllchem1 + llz1; 
      vsllchem1 = vsllchem1 + 2 * s ^ 2 * (llzero - llz1);
    } 

    if (TN[n] >= 0) {
      vnTN = vnTN + 1;
      s = chem2ae + chem2re*bach[n,5];
      z2 = (TN[n] - bach[n,5]) / s;
      vsumz2 = vsumz2 + z2;
      vsumz22 = vsumz22 + z2 * z2;
      if (!is_nan(lastz2)) {
        vsumz2a = vsumz2a + z2 * lastz2;
        vsumz2b = vsumz2b + z2 + lastz2;
      }
      lastz2 = z2;
      llz2 = student_t_lpdf(z2 | 7, 0, 1);
      vllchem2 = vllchem2 + llz2; 
      vsllchem2 = vsllchem2 + 2 * s ^ 2 * (llzero - llz2); 
    }

  } // validation calcs

  lastz1 = 0; // illegal to return nan
  lastz2 = 0;

  vmeanz1 = vnTP>2 ? vsumz1/vnTP : 0;
  vmeanz2 = vnTN>2 ? vsumz2/vnTN : 0;
  vauto1 = vnTP>2 ? (vsumz1a-vmeanz1*vsumz1b+(vnTP-1)*vmeanz1^2)/(vsumz12-vnTP*vmeanz1^2) : 0;
  vauto2 = vnTN>2 ? (vsumz2a-vmeanz2*vsumz2b+(vnTN-1)*vmeanz2^2)/(vsumz22-vnTN*vmeanz2^2) : 0;

  vllchem1 = -vllchem1;
  vllchem2 = -vllchem2;
  vlltotal = vllchem1 + vllchem2;
  vgrmse1 = vnTP>2 ? sqrt(vsllchem1/vnTP) : 0;
  vgrmse2 = vnTN>2 ? sqrt(vsllchem2/vnTN) : 0;

}

