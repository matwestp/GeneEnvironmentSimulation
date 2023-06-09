//Simulation für das Gene Paper: Interaktion zwischen beobachteter und unbeobachteter Heterogenität

cd "C:\Users\Win7ADM\Documents\GitHub\GeneEnvironmentSimulation"
graph set window fontface "Century"
set scheme plotplainblind
*******************************************************************************

	local varU =4
	clear 
	set obs 100000

	**# Generate Variables
	// Unobserved heterogeneity
	loc corr =.4
	mat C = 1,`corr' \ `corr',1
	drawnorm eps0 eps1, corr(C)

	//Endowment (Genes)
// 	gen 	G =eps0>0 
gen G =rnormal()>0								// Leute mit positiven (unbeobachteten) Y-Werten haben gute Gene corr(eps0,eps1)>0
	la var 	G "Advantageous genetic environment"

	
	//Potential Outcomes 
	local ATE_D_G0				=.3
	local ATE_D_interaction 	=1.5
	local ATE_G_D0				=.5
	local intercept				=0 
	
	gen Y1 =(`intercept' +`ATE_D_G0') 	+(`ATE_G_D0'+`ATE_D_interaction')*G + eps1
	gen Y0 =`intercept' + 				  `ATE_G_D0'*G 						+ eps0

	//Instrument
	gen Z =rnormal()>0

	//Environment (Education); the Treatment
	gen D0 	=(Y1-Y0)-.53+1*G >0
	gen D1 	=(Y1-Y0)-.53+1*G >(`varU'*(1-.75*G))
	gen D 	=D0 +Z*(D1-D0)

	//Define Compliers
	gen comp =(D0==1)*(D1==0)


	bys G: su eps0 if D==1

	gen ITE =Y1-Y0 
	su 	ITE if D==1 & comp==1
	loc LATE_true = `r(mean)'

	bys G: su ITE if D==1

	bys G: su comp 

	//Observation rule 
	gen Y =Y0 +D*(Y1-Y0)
	// OLS is biased because cov(eps0,D)!=0
	// IV is unbiased because cov(eps0,Z)==0
	*********************************************************************

	**# Analyzing Variables 


	qui{ // Theorem IV weights (Loken, Mogstad, Wiswall, 2012 AEJ:Applied) --- Answers the question of what weight the sub-LATEs in G=1 and G=0 receive for the overall LATE
	ivregress 2sls Y (D=Z)  if G==0
	loc LATE_G0 =_b[D]
	ivregress 2sls Y (D=Z) if G==1
	loc LATE_G1 =_b[D]

	ivregress 2sls Y (D=Z) if G==1
	loc Z1 =_b[D]
		cap drop DG?
	gen DG1 =D*G
	ivregress 2sls Y (D=Z) if G==0
	loc Z0 =_b[D]
	
	gen DG0 =D*(G==0)
	ivregress 2sls DG0 (D=Z)
	loc weightG0 =_b[D]

	loc weightedLATE `Z1'*(1-`weightG0')+`Z0'*(`weightG0')
	
	ivregress 2sls Y (D=Z) 
	loc LATE =_b[D]
	noi di "weighted LATE: " _col(20) `weightedLATE'
	noi di "LATE:" _col(20) `LATE'
	noi di "True LATE:" _col(20) `LATE_true'
	noi di _n(1) "Note: in this LATE, two mechanisms operate:"
	noi di _col(5) "1. G affects who are the compliers are (measured in terms of U_D, referred to as selection effect)." 
	noi di _col(5) "2. G affecs the outcomes differently (for the same individuals, i.e. conditional on U{sub:D})"
	noi di _n(1) _col(3) "==> Now we will decompose both effects..."
	}


	qui{
	gen V =ITE // In the case of control variables, one need to fix them as U_D varies 
	la var V "unobserved gains"
	gen U_D =. 
	foreach g of numlist 0 1{
		su V if G==`g'
		replace U_D =normal(-(V-`r(mean)')/`r(sd)') if G==`g'
		su ITE if D==1 &G==`g' &comp==1
		glo mATT`g' =`r(mean)'
		su comp if G==`g'
		glo comp`g' =`r(mean)'
	}
	}
	
	**# Efefct decomposition

	gen eval =_n/100 if _n<=100 
	
	qui{ // Estimate the weights. Theorem IV weights (Loken, Mogstad, Wiswall, 2012 AEJ:Applied) --- Answers the question of what weight the sub-LATEs in each U_D bin receive for the overall LATE by G 
			cap drop wG? 
		foreach g of numlist 0 1{
			gen wG`g' =. 
			foreach num of numlist .05(.05).95{
				di `num'
					cap drop DUD`=floor(`num'*100)'
				gen DUD`=floor(`num'*100)' = D*inrange(U_D,`=`num'-.025',`=`num'+.025')

				cap ivregress 2sls DUD`=floor(`num'*100)' (D=Z) if G==`g'
				di "`num': "_b[D]
				if  _rc==0 	replace wG`g' =_b[D] 	if inrange(eval,`num'-.001,`num'+.001)
				else 		replace wG`g' =0 		if inrange(eval,`num'-.001,`num'+.001)
			}
		}
	}	

//	Estimate/compute the true MTEs 
	lpoly V U_D if G==1, gen(MTE1) at(eval) nogr
	lpoly V U_D if G==0, gen(MTE0) at(eval) nogr

	la var MTE1 "MTE, G=1" 
	la var MTE0 "MTE, G=0"
    
//	Estimate the complier weights 
	lpoly comp U_D if G==0, gen(wLATE0) at(eval) nogr
	lpoly comp U_D if G==1, gen(wLATE1) at(eval) nogr

	la var wLATE0 "{&omega}{sub:LATE}, G=0"
	la var wLATE1 "{&omega}{sub:LATE}, G=1"

//	the lines for the LATE
	gen LATE1 =$mATT1 if !mi(eval) &wLATE1>0
	gen LATE0 =$mATT0 if !mi(eval) &wLATE0>0
	
	la var LATE1 "LATE, G=1"
	la var LATE0 "LATE, G=0"
	
	la var eval "U{sub:E}"
	
	qui{ // Determine points at which to display the effects 
	su eval if wLATE1>0
	loc evalmin1 =`r(min)'
	loc evalmax1 =`r(max)'
	su eval if wLATE0>0
	loc evalmin0 =`r(min)'
	loc evalmax0 =`r(max)'
	local evalmean =floor((min(`evalmin0',`evalmin1')+max(`evalmax0',`evalmax1'))/2*100)
	di `evalmean'
	
	if abs(50-`evalmean')<20	loc ev2 =10
	else 						loc ev2 =50
	}
	
	local TrueGEtext 	"True G x E interaction"
	local IVGEtext 		"IV estimate of G x E interaction"
	
	gen 	MTElab ="`TrueGEtext'" 	if _n==`ev2'
	replace MTElab ="`IVGEtext'" 	if _n==`evalmean'
	
	gen mean =(MTE1+MTE0)/2 if _n==`ev2'
	su LATE1
	loc LATE1 =`r(mean)'
	su LATE0
	loc LATE0 =`r(mean)'
	di "inlist(_n,`evalmean',`ev2')"
	
	gen eval1 =eval if inlist(_n,`evalmean',`ev2')
	replace mean =(`LATE1'+`LATE0')/2 if _n==`evalmean'
	gen effect1 =MTE1 if _n==`ev2'
	replace effect1 =`LATE1' if _n==`evalmean'
	gen effect0 =MTE0 if _n==`ev2'
	replace effect0 =`LATE0' if _n==`evalmean'
	colorpalette plottig
//  ananas sky turquoise reddish vermillion sea orangebrown
// 	tw (area wG0 eval, col("97 156 255 %30")) (area wLATE0 eval, yaxis(2) col("0 176 246 %30"))
// 	tw (li wG1 eval) (li wLATE1 eval, yaxis(2))
	global col1 "97 156 255"
// 	global col2 "0 176 246"
// 	global col1 "0 188 216"
	global col2 "0 192 175"
//	MTE plot: 
//Presentation, slide 1
	tw (li MTE1 MTE0 eval, lc("$col1" none) lw(.6 =)) (area wLATE1 wLATE0 eval, lp(dash =) yaxis(2) col(none =)) (li LATE1 LATE0 LATE1 LATE0 eval, lw(1 = = =) lp(solid = dot =) lc(none = = = )) (rcap effect? eval1, mlab(MTElab) col(none =)) (sc mean eval1, mlab(MTElab) msize(0)  mlabcol(none =)), plotr(lc(none)) legend(pos(6) order(1 3 9 2 4 10) cols(3)) ytitle("`Pr(Complier)'", axis(2)) ytitle("Effect") name(Gr`=`varU'*10', replace) xtitle("`:var label eval'") /*title("{&sigma}{sub:U}=`varU'")*/
	gr export "Setting_1.pdf", replace 
//	slide 2
	tw (li MTE1 MTE0 eval, lc("$col1" "$col2") lw(.6 =)) (area wLATE1 wLATE0 eval, lp(dash =) yaxis(2) col(none =)) (li LATE1 LATE0 LATE1 LATE0 eval, lw(1 = = =) lp(solid = dot =) lc(none = = = )) (rcap effect? eval1, mlab(MTElab) col(none =)) (sc mean eval1, mlab(MTElab) msize(0)  mlabcol(none =)), plotr(lc(none)) legend(pos(6) order(1 3 9 2 4 10) cols(3)) ytitle("`Pr(Complier)'", axis(2)) ytitle("Effect") name(Gr`=`varU'*10', replace) xtitle("`:var label eval'") /*title("{&sigma}{sub:U}=`varU'")*/
	gr export "Setting_2.pdf", replace 
//	slide 3
	tw (li MTE1 MTE0 eval, lc("$col1" "$col2") lw(.6 =)) (area wLATE1 wLATE0 eval, lp(dash =) yaxis(2) col(none =)) (li LATE1 LATE0 LATE1 LATE0 eval, lw(1 = = =) lp(solid = dot =) lc(none = = = )) (rcap effect0 effect1 eval1 if _n==10, col(black) mlab(MTElab)) (sc mean eval1 if _n==10, mlab(MTElab) msize(0)  mlabcol(black none)), plotr(lc(none)) legend(pos(6) order(1 3 9 2 4 10) cols(3)) ytitle("`Pr(Complier)'", axis(2)) ytitle("Effect") name(Gr`=`varU'*10', replace) xtitle("`:var label eval'") /*title("{&sigma}{sub:U}=`varU'")*/
	gr export "Setting_3.pdf", replace 
//	slide 4 
	tw (li MTE1 MTE0 eval, lc("$col1" "$col2") lw(.6 =)) (area wLATE1 wLATE0 eval, lp(dash =) yaxis(2) col("$col1 %50" none)) (li LATE1 LATE0 LATE1 LATE0 eval, lw(1 = = =) lp(solid = dot =) lc("$col1 %50"  none black none)) (rcap effect0 effect1 eval1 if _n==10, col(black) mlab(MTElab)) (sc mean eval1 if _n==10, mlab(MTElab) msize(0)  mlabcol(black none)), plotr(lc(none)) legend(pos(6) order(1 3 9 2 4 10) cols(3)) ytitle("`Pr(Complier)'", axis(2)) ytitle("Effect") name(Gr`=`varU'*10', replace) xtitle("`:var label eval'") /*title("{&sigma}{sub:U}=`varU'")*/
	gr export "Setting_4.pdf", replace 
//	slide 5 
	tw (li MTE1 MTE0 eval, lc("$col1" "$col2") lw(.6 =)) (area wLATE1 wLATE0 eval, lp(dash =) yaxis(2) col("$col1 %50" "$col2 %50")) (li LATE1 LATE0 LATE1 LATE0 eval, lw(1 = = =) lp(solid = dot =) lc("$col1 %50"  "$col2 %50" black =)) (rcap effect0 effect1 eval1 if _n==10, col(black) mlab(MTElab)) (sc mean eval1 if _n==10, mlab(MTElab) msize(0)  mlabcol(black none)), plotr(lc(none)) legend(pos(6) order(1 3 9 2 4 10) cols(3)) ytitle("`Pr(Complier)'", axis(2)) ytitle("Effect") name(Gr`=`varU'*10', replace) xtitle("`:var label eval'") /*title("{&sigma}{sub:U}=`varU'")*/
	gr export "Setting_4.pdf", replace 
	
	tw (li MTE1 MTE0 eval, lc("$col1" "$col2") lw(.6 =)) (area wLATE1 wLATE0 eval, lp(dash =) yaxis(2) col("$col1 %50" "$col2 %50")) (li LATE1 LATE0 LATE1 LATE0 eval, lw(1 = = =) lp(solid = dot =) lc("$col1" "$col2")) (rcap effect? eval1, mlab(MTElab)) (sc mean eval1, mlab(MTElab) msize(0)), plotr(lc(none)) legend(pos(6) order(1 3 9 2 4 10) cols(3)) ytitle("`Pr(Complier)'", axis(2)) ytitle("Effect") name(Gr`=`varU'*10', replace) xtitle("`:var label eval'") /*title("{&sigma}{sub:U}=`varU'")*/
	gr export "Setting_final.pdf", replace 
	
	
	qui{ // Bar Plot to visualize effect sizes
		gen deffect =effect1-effect0
		la var deffect "Effect size"
		gen mideffect =mi(deffect)
		bys mideffect: gen bareval =_n if mideffect==0
		sort eval
		su bareval if MTElab=="`TrueGEtext'"
		loc evalOutEff =`r(mean)'
		su bareval if MTElab=="`IVGEtext'"
		loc evalIVEff =`r(mean)'
		su deffect
		tw (bar deffect bareval, barw(.75) col(black)), ylabel(`=ceil(`r(min)')'(.5)`=ceil(`r(max)')') yline(0) xlabel(`evalOutEff' "`TrueGEtext'" `evalIVEff' "`IVGEtext'") xtitle("") plotr(lc(none)) 
		gr di, xsize(2.5) name(bar`=`varU'*10', replace)
		gr export "Effect_comparison_`=`varU'*10'.pdf", replace 
	}
	gen dMTE =MTE1 - MTE0 
	su dMTE
	glo outcome_effect =`r(mean)'
	
	noi di "1. Outcome effect: " _col(21)$outcome_effect
	noi di "2. Selection effect: "_col(21)`LATE' -$outcome_effect
	noi di _col(5) "(Difference in complying probability is: `=round(abs($comp1-$comp0),.0001)')"
// 	di  $mATT1 - $mATT0
	noi di _n(1) "G x E model is:"
	noi di _col(5) "ivregress 2sls Y (D 1.D#1.G=Z 1.Z#1.G) 1.G"
	noi di _n(1) "`IVGEtext' is" 
	noi di _col(5) "_b[1.D#1.G]"

// exit

gen ZG1 =Z*G
ivregress 2sls Y (D 1.D#1.G=Z 1.Z#1.G) 1.G
grc1leg Gr40 bar40, name(combined, replace) 
gr di combined, ysize(3) xsize(6) scale(1.2)
gr export "Simulation_results.pdf", replace 

********************************************************************************

*synchronize file 
file close _all
file open git using mygit.bat, write replace 
file write git "git remote add origin " `"""' "https://github.com/matwestp/GeneEnvironmentSimulation.git" `"""' _n
 file write git "git add --all" _n
 file write git "git commit -m "
 file write git `"""' "change" `"""' _n 			// choose name for change 
 file write git "git push" _n
file close git

! mygit.bat