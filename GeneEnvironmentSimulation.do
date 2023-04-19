//Simulation für das Gene Paper: Interaktion zwischen beobachteter und unbeobachteter Heterogenität

cd "C:\Users\Win7ADM\Documents\GitHub\GeneEnvironmentSimulation"

*******************************************************************************
clear 
set obs 10000

mat C = 1,-.0\-.0,1

drawnorm eps0 eps1, corr(C)

gen Y1 =.2*eps1
gen Y0 =.2*eps0

gen Z =rnormal()>0

//Environment (Education); the Treatment
gen D0 =(Y1-Y0)>0
gen D1 =(Y1-Y0)>1
gen D =D0 +Z*(D1-D0)

gen comp =(D0==1)*(D1==0)
gen ITE =Y1-Y0 
su ITE if D==1
//Endowment (Genes)
gen G =eps0>0

bys G: su ITE if D==1

bys G: su comp 


gen Y =Y0 +D*(Y1-Y0)

reg Y D 

di "aUtO"

qui{
reg Y 1.D#i.G
loc ATT_G0 =_b[1.D#0.G]
loc ATT_G1 =_b[1.D#1.G]

su G if D==1 
loc Gmean =`r(mean)'

reg Y D 
loc ATT =_b[D]
noi di "Gewichter ATT: " _col(20) `ATT_G0'*(1-`Gmean')+`ATT_G1'*`Gmean'
noi di "ATT:" _col(20) `ATT'
}


gen V =ITE 
la var V "unobserved gains"
gen U_D =. 
foreach g of numlist 0 1{
	su V if G==`g'
	replace U_D =normal(-(V-`r(mean)')/`r(sd)') if G==`g'
	su ITE if D==1 &G==`g'
	loc mATT`g' =`r(mean)'
}

tw (lpoly V U_D if G==1, lc(blue)) (lpoly V U_D if G==0, lc(red)) (lpoly comp U_D if G==1, lp(dash) yaxis(2) lc(blue)) (lpoly comp U_D if G==0, lp(dash) yaxis(2) lc(red)), yline(`mATT1' `mATT0', lc(blue red))

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
