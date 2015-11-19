# First simulate the marker effects
using Distributions

function simDat(nObs,nLoci,bMean,bStd,resStd)
X=[ones(nObs,1) sample([0,1,2],(nObs,nLoci))]
b=rand(Normal(bMean,bStd),size(X,2))
y=X*b+rand(Normal(0.0,resStd),nObs)
return(y,X,b)
end

nObs =10;
nLoci =25;
bMean =0.0;
bStd =0.5;
resStd =1.0;
res=simDat(nObs,nLoci,bMean,bStd,resStd);
y=res[1];
X=res[2];
b=res[3];

# Convert the data objects to my notation
M=X[:,2:end]  
X=X[:,1]
Im=eye(nLoci)
sigmaSqE=1
sigmaSqA=.2
lambdaSNP=sigmaSqE/sigmaSqA
M=[M;ones(nLoci)']
X=[X;1]
nObs+=1
y=[X M]*b+rand(Normal(0.0,resStd),nObs)
oldM=M
#M=M-1
#M=2-M

# Fit the marker effects model (SNP-BLUP)
SNPmmeLhs=[X'X X'M
           M'X (M'M+Im*lambdaSNP)]
SNPmmeRhs=[X'y
           M'y]
SNPinvMME=inv(SNPmmeLhs)
SNPMarkersoln=SNPmmeLhs\SNPmmeRhs


SNP_EBVsoln=M*SNPMarkersoln[2:end]
SNPC22=SNPinvMME[2:end,2:end]
varuhat=M*(Im*sigmaSqA - SNPC22*sigmaSqE)*M'
varu=M*M'*sigmaSqA
pev_uhat=M*(SNPC22*sigmaSqE)*M'

relSNP_EBV=diag(varuhat)./diag(varu)




# Fit the breeding value model (GBLUP)
G=varu
Z=eye(nObs)
Ginv=inv(G)

mmeLhs=[X'X X'Z
        Z'X (Z'Z+Ginv*sigmaSqE)]
mmeRhs=[X'y
        Z'y]
invMME=inv(mmeLhs)
EBVsoln=mmeLhs\mmeRhs
C22=invMME[2:end,2:end]
rel_EBV=diag(G-C22)./diag(G)

println("Original Parameterizations")
[EBVsoln[2:end] rel_EBV SNP_EBVsoln relSNP_EBV]


M=oldM-1

# Fit the marker effects model (SNP-BLUP)
SNPmmeLhs=[X'X X'M
           M'X (M'M+Im*lambdaSNP)]
SNPmmeRhs=[X'y
           M'y]
SNPinvMME=inv(SNPmmeLhs)
SNPMarkersoln=SNPmmeLhs\SNPmmeRhs


SNP_EBVsoln=M*SNPMarkersoln[2:end]
SNPC22=SNPinvMME[2:end,2:end]
varuhat=M*(Im*sigmaSqA - SNPC22*sigmaSqE)*M'
varu=M*M'*sigmaSqA
pev_uhat=M*(SNPC22*sigmaSqE)*M'

relSNP_EBV=diag(varuhat)./diag(varu)




# Fit the breeding value model (GBLUP)
G=varu
Z=eye(nObs)
Ginv=inv(G)

mmeLhs=[X'X X'Z
        Z'X (Z'Z+Ginv*sigmaSqE)]
mmeRhs=[X'y
        Z'y]
invMME=inv(mmeLhs)
EBVsoln=mmeLhs\mmeRhs
C22=invMME[2:end,2:end]
rel_EBV=diag(G-C22)./diag(G)

println("M=M-1 Parameterizations")
[EBVsoln[2:end] rel_EBV SNP_EBVsoln relSNP_EBV]


M=2-oldM
 
 
# Fit the marker effects model (SNP-BLUP)
SNPmmeLhs=[X'X X'M
           M'X (M'M+Im*lambdaSNP)]
SNPmmeRhs=[X'y
           M'y]
SNPinvMME=inv(SNPmmeLhs)
SNPMarkersoln=SNPmmeLhs\SNPmmeRhs


SNP_EBVsoln=M*SNPMarkersoln[2:end]
SNPC22=SNPinvMME[2:end,2:end]
varuhat=M*(Im*sigmaSqA - SNPC22*sigmaSqE)*M'
varu=M*M'*sigmaSqA
pev_uhat=M*(SNPC22*sigmaSqE)*M'

relSNP_EBV=diag(varuhat)./diag(varu)




# Fit the breeding value model (GBLUP)
G=varu
Z=eye(nObs)
Ginv=inv(G)

mmeLhs=[X'X X'Z
        Z'X (Z'Z+Ginv*sigmaSqE)]
mmeRhs=[X'y
        Z'y]
invMME=inv(mmeLhs)
EBVsoln=mmeLhs\mmeRhs
C22=invMME[2:end,2:end]
rel_EBV=diag(G-C22)./diag(G)

println("M=2-M Parameterizations")
[EBVsoln[2:end] rel_EBV SNP_EBVsoln relSNP_EBV]
 