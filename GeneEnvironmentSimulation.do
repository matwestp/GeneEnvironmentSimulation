//Simulation für das Gene Paper: Interaktion zwischen beobachteter und unbeobachteter Heterogenität

cd "C:\Users\Win7ADM\Documents\GitHub\GeneEnvironmentSimulation"

*******************************************************************************
clear 
set obs 10000

**# Generate Variables
// Unobserved heterogeneity
loc corr =-.2
mat C = 1,`corr' \ `corr',1
drawnorm eps0 eps1, corr(C)

//Potential Outcomes 
gen Y1 =.2*eps1
gen Y0 =.2*eps0

//Instrument
gen Z =rnormal()>0

//Environment (Education); the Treatment
gen D0 =(Y1-Y0)>0
gen D1 =(Y1-Y0)>1
gen D =D0 +Z*(D1-D0)

//Define Compliers
gen comp =(D0==1)*(D1==0)

//Endowment (Genes)
gen G =eps0>0

gen ITE =Y1-Y0 
su ITE if D==1
loc ATT_true = `r(mean)'

bys G: su ITE if D==1

bys G: su comp 

//Observation rule 
gen Y =Y0 +D*(Y1-Y0)
// OLS is biased because cov(eps0,D)!=0
// IV is unbiased because cov(eps0,Z)==0
*********************************************************************

**# Analyzing Variables 


qui{
*reg Y 1.D#i.G
// loc ATT_G0 =_b[1.D#0.G]
// loc ATT_G1 =_b[1.D#1.G]
ivregress 2sls Y (D=Z)  if G==0
loc ATT_G0 =_b[D]
ivregress 2sls Y (D=Z) if G==1
loc ATT_G1 =_b[D]

su G if D==1 
loc Gmean =`r(mean)'

reg Y D 
ivregress 2sls Y (D=Z) 
loc ATT =_b[D]
noi di "Gewichter ATT: " _col(20) `ATT_G0'*(1-`Gmean')+`ATT_G1'*`Gmean'
noi di "ATT:" _col(20) `ATT'
noi di "True ATT:" _col(20) `ATT_true'
noi di "Note: in this ATT, two mechanisms operate: 1. effects on complying probability. 2. Effect on outcomes conditional on complying probability"
}

exit

gen V =ITE 
la var V "unobserved gains"
gen U_D =. 
foreach g of numlist 0 1{
	su V if G==`g'
	replace U_D =normal(-(V-`r(mean)')/`r(sd)') if G==`g'
	su ITE if D==1 &G==`g'
	glo mATT`g' =`r(mean)'
	su comp if G==`g'
	glo comp`g' =`r(mean)'
}



gen eval =_n/100 if _n<=100 

lpoly V U_D if G==1, gen(MTE1) at(eval)
lpoly V U_D if G==0, gen(MTE0) at(eval)

la var MTE1 "MTE, G=1" 
la var MTE0 "MTE, G=0"

gen LATE1 =$mATT1 if !mi(eval) &eval<=$comp1
gen LATE0 =$mATT0 if !mi(eval) &eval<=$comp0

la var LATE1 "LATE, G=1"
la var LATE0 "LATE, G=0"

lpoly comp U_D if G==0, gen(wLATE0) at(eval)
lpoly comp U_D if G==1, gen(wLATE1) at(eval)

la var wLATE0 "{&omega}{sub:LATE}, G=0"
la var wLATE1 "{&omega}{sub:LATE}, G=1"

la var eval "U{sub:D}"

tw (li MTE1 MTE0 eval, lc(blue red) lw(.6 =)) (li wLATE1 wLATE0 eval, lp(dash =) yaxis(2) lc(blue red)) (li LATE1 LATE0 eval, lw(1 =) lc(blue red)) , plotr(lc(none)) legend() ytitle("Weight", axis(2)) ytitle("Effect")

gen dMTE =MTE1 - MTE0 
su dMTE

di  $mATT1 - $mATT0
bys G: ivregress 2sls Y (D=Z) 


reg Y D 

su V if D==1

********************************************************************************

*synchronize file (cd muss immer lokal auf den github ordner eingestellt sein)
file close _all
file open git using mygit.bat, write replace 
file write git "git remote add origin " `"""' "hhttps://github.com/matwestp/GeneEnvironmentSimulation.git" `"""' _n
file write git "git add --all" _n
file write git "git commit -m "
file write git `"""' "change" `"""' _n 			// choose name for change 
file write git "git push" _n
file close git

! mygit.bat
