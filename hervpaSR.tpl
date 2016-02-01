//><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
//Programer: UBC ADMB May Workshop													 
//Date:	May 12-16, 2013														 
//Purpose:Herring VPA.											 
//Notes: 				 
//							 
//																 
//><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>

DATA_SECTION
	init_int syr;
	init_int eyr;
	init_int nage;
	init_int agefv;
	init_vector surv(1,nage);
	init_vector mat(1,nage);
	init_vector wa(1,nage);
	init_vector yt(syr,eyr);
	init_matrix cat(syr,eyr,1,nage);
	init_int eof;
	int iter;
	!!iter=0;
	LOCAL_CALCS
		if(eof!=999)
		{
			cout<<"Error reading data.\n Fix it."<<endl;
			ad_exit(1);
		}
	END_CALCS
	vector age(1,nage);
	vector fa(1,nage);
	matrix vat(syr,eyr,1,nage);
	LOCAL_CALCS
		age.fill_seqadd(1,1);
		fa=elem_prod(wa,mat);
	END_CALCS
PARAMETER_SECTION
	init_bounded_number U(0,1);
	init_number log_ptau;
	!!log_ptau=1./log(0.3*0.3);
	init_number log_psig(2);
	!!log_psig=1./log(0.3*0.3);
	init_number log_ro(2);
	!!log_ro=log(1200);
	init_number log_cr(2);
	!!log_cr=log(6);
	objective_function_value nll;
	
	number sig;
	number tau;
	number q;
	number ro;
	number cr;
	vector selt(1,nage);
	vector rt_resid(syr,eyr-1);
	vector yt_resid(syr,eyr);
	vector SBt(syr,eyr);
	vector Ut(syr,eyr);
	vector Rt(syr,eyr);
	matrix Nt(syr,eyr,1,nage);

PROCEDURE_SECTION
	initialization();
	statedynamics();
	observation_model();
	stock_recruit_model();
	objective_function();

	if(mceval_phase())
	{ 
		forecast();
		mcmc_output();
	}
	if(last_phase())
	{
		forecast();
	} 

FUNCTION initialization
	Nt(eyr)(agefv,nage)=cat(eyr)(agefv,nage)/U;
	selt(agefv,nage)=1.;
	Ut(eyr)=U;

FUNCTION statedynamics
	// do the stitching for fully vulnerable individuals and
	//1 age younger than fully vulnerable
	for(int ti=eyr-1;ti>=syr;ti--)
	{
		for(int a=agefv-1;a<=nage-1;a++)
		{
			Nt(ti,a)=Nt(ti+1,a+1)/surv(a)+cat(ti,a);
		}
		Ut(ti)=sum(cat(ti)(agefv,nage-1))/sum(Nt(ti)(agefv,nage-1));
		Nt(ti,nage)=cat(ti,nage)/Ut(ti);
	}

	//stitch for younger ages getting vat from a 
	// 4 year window
	int twindow=4;
	for(int a=agefv-1;a>=2;a--)
	{
		dvariable va=0;
		int nva=0;
		for(int ti=eyr-twindow;ti<=eyr-1;ti++)
		{
			va+=(cat(ti,a)/Nt(ti,a))/Ut(ti);
			nva++;
		}

		selt(a)=va/nva;
		selt(a)+=1.e-8;
		Nt(eyr,a)=cat(eyr,a)/(U*selt(a));
		for(int ti=eyr;ti>=syr+1;ti--)
		{
			Nt(ti-1,a-1)=Nt(ti,a)/surv(a-1)+cat(ti-1,a-1);
		}
	}
	//calculate the vulnerability for the first age and
	//the numbers in the terminal year in the first age
	dvariable va=0;
	int nva=0;
	for(int ti=eyr-twindow;ti<=eyr-1;ti++)
	{
		va+=(cat(ti,1)/Nt(ti,1))/Ut(ti);
		nva++;
	}
	selt(1)=va/nva;
	Nt(eyr,1)=cat(eyr,1)/(U*selt(1));
	vat=1;
	for(int ti=syr;ti<=eyr;ti++)
	{
		for(int a=1;a<=nage;a++)
		{
			if(Nt(ti,a)>0)vat(ti,a)=(cat(ti,a)/value(Nt(ti,a)))/value(Ut(ti));
		}
	}
	Rt=column(Nt,1);

FUNCTION observation_model
	SBt=Nt*fa;
	yt_resid=log(yt)-log(SBt);
	q=mfexp(mean(yt_resid));
	yt_resid-=mean(yt_resid);

FUNCTION stock_recruit_model
	ro=mfexp(log_ro);
	cr=mfexp(log_cr)+1.;
	dvar_vector lxo(1,nage);
	lxo(1)=1;
	for(int i=2;i<=nage;i++)lxo(i)=lxo(i-1)*surv(i-1);
	lxo(nage)/=(1.-surv(nage));
	dvariable phieo=lxo*fa;
	dvariable so=cr/phieo;
	dvariable beta=(cr-1.)/(ro*phieo);
	dvar_vector sbt=Nt*fa;
	dvar_vector nmr=so*sbt(syr,eyr-1);
	dvar_vector den=(1.+beta*sbt(syr,eyr-1));
	dvar_vector tmp_rt=elem_div(nmr,den);
	dvar_vector rt=--column(Nt.sub(syr+1,eyr),1);
	rt_resid=log(tmp_rt)-log(rt);
	//cout<<rt_resid<<endl;
	//ad_exit(1);
FUNCTION objective_function 
	tau=sqrt(1./mfexp(log_ptau));
	sig=sqrt(1./mfexp(log_psig));
	
	nll=dnorm(yt_resid,tau)+dnorm(rt_resid,sig);

FUNCTION mcmc_output
	if(iter==0)
	{
		//ofstream ofs("refpar.mcmc");
		//ofs<<"fmsy\t bmsy\t msy\t b/bmsy\t f/fmsy"<<endl;
		//ofs<<"r\t k\t q\t sig\t"<<endl;
	}
	iter++;
	//double fratio=value(ft(eyr)/fmsy);
	//double bratio=value(Nt(eyr)*wa/bmsy);
	//ofstream ofs("refpar.mcmc",ios::app);
	//ofs<<fmsy<<"\t"<<bmsy<<"\t"<<msy<<"\t"<<bratio<<"\t"<<fratio<<endl;

FUNCTION forecast


REPORT_SECTION
	REPORT(U);	
	REPORT(tau);	
	REPORT(sig);	
	REPORT(cr);	
	REPORT(ro);	
	REPORT(q);
	REPORT(Ut);	
	REPORT(SBt);	
	REPORT(Rt);
	REPORT(yt_resid);
	REPORT(Nt);	
	REPORT(vat);	
	REPORT(cat);	

TOP_OF_MAIN_SECTION
	
	time(&start);
	arrmblsize = 50000000;
	gradient_structure::set_GRADSTACK_BUFFER_SIZE(1.e7);
	gradient_structure::set_CMPDIF_BUFFER_SIZE(1.e7);
	gradient_structure::set_MAX_NVAR_OFFSET(5000);
	gradient_structure::set_NUM_DEPENDENT_VARIABLES(5000);


GLOBALS_SECTION
	/**
	\def REPORT(object)
	Prints name and value of \a object on ADMB report %ofstream file.
	*/
	#undef REPORT
	#define REPORT(object) report << #object "\n" << object << endl;

	#include <admodel.h>
	#include <time.h>
	#include <contrib.h>//IF you have ADMB-11
	//#include<stats.cxx>//If you have ADMB-10 and make sure stats.cxx is in your working directory
	time_t start,finish;
	long hour,minute,second;
	double elapsed_time;

FINAL_SECTION
	time(&finish);
	elapsed_time=difftime(finish,start);
	hour=long(elapsed_time)/3600;
	minute=long(elapsed_time)%3600/60;
	second=(long(elapsed_time)%3600)%60;
	cout<<"*******************************************"<<endl;
	cout<<"--Start time: "<<ctime(&start)<<endl;
	cout<<"--Finish time: "<<ctime(&finish)<<endl;
	cout<<"--Runtime: ";
	cout<<hour<<" hours, "<<minute<<" minutes, "<<second<<" seconds"<<endl;
	cout<<"*******************************************"<<endl;


