#include "DASpec.h"
// #include "cpgplot.h"
#include "fstream"

int main(int argc, char const *argv[])
{
    int num = 7415;
    double x[num], y[num];
    double x1[num], y1[num], err1[num];
    double x2[num], y2[num];
    float plotx[num], ploty[num];
    
    std::ifstream file("fetemplate_no3");
    for (int i = 0; i < num; ++i)
    {
        file >> x[i] >> y[i];
        x1[i] = x[i];
        y1[i] = 0.0;
    //    x2[i] = x[i];
    //    y2[i] = 0.0;
        err1[i] = 1.0e5;
    //    //std::cout << x[i] << " " << y[i] << std::endl;
    }
    
    //balmer_continuum balmer = balmer_continuum();
    //double p[] = {1.0, 0.5};
    //balmer.calc(num, x, y, p);
    //
    powerlaw a = powerlaw(5100);
    
    //template_spec a = template_spec(num, x, y, "gaussian");
    //double p[] = {20.0, 2000.0, 1000.0};
    double p[] = {20.0, -2.2};
    a.calc(num, x1, y1, p);
    a.info();
    
    //template_spec_reddened b = template_spec_reddened(num, x, y, "gaussian");
    //double p1[] = {2.0, 100.0, 1000.0, 0.5};
    //b.calc(num, x2, y2, p1);
    //b.info();
    
    //template_spec
    
    //for (int i = 0; i< num; ++i)
    //{
    //    x[i] = 4500.0 + i;
    //    y[i] = 0.0;
    //    err[i] = 0.01;
    //}
        
	//component *func;
    compcontainer model;
    //model.add(new template_spec(num, x, y, "gaussian"));
    model.add(new powerlaw(5001));
    //model.add(new line_gaussian(4861));
    //model.add(new line_gaussian(5007));
    
    //template_spec a = template_spec(num, x, y, "gh4");
    //a.info();
    
    //std::exit(-1);
    //model.info();
    //model.addfix(1, 1, 0.3);
    //model.addfix(1, 2, 2.0);
    //model.addtie(0, 1, 1, 1, "ratio", 2.0);
    //model.addtie(0, 2, 1, 2, "offset", 10.0);
    //model.addtie_profile(3, 2);
    //double p0[3] = {50.0, 1000.0, 0.0};
    double p0[2] = {1.0, -1.0};
    //model.calc(num, x, y, p);
    //model.info();
    //double ptot[model.npar];
    //model.pars2l(p, ptot);
    //std::cout << "set parameters to:" << std::endl;
    //for (int i = 0; i < model.npar; ++i) std::cout << ptot[i] << " ";
    //std::cout << std::endl;
    
    //compcontainer model1;
    //model1.add(new powerlaw(5001));
    //model1.add(new line_gaussian(4861));
    //model1.add(new line_gaussian(5007));
    //model1.addfix(1, 1, 0.3);
    //model1.addtie_profile(3, 2);
    
    curvefit fit;
    //double p0[5] = {-1.0, 10.0, 6000.0, 0.0, 0.0};
    fit.setdata(num, x1, y1, err1);
    fit.setinitp(p0);
    fit.setmodel(&model);
    fit.setlimit(1, 0.0, 0);
    fit.setlimit(1, 100.0, 1);
    fit.setlimit(2, -5.0, 0);
    fit.setlimit(2, 5.0, 1);
    //fit.setlimit(1, 0.0, 0);
    //fit.setlimit(1, 1.0e3, 1);
    //fit.setlimit(2, 100.0, 0);
    //fit.setlimit(2, 5.0e3, 1);
    //fit.setlimit(3, -1.0e3, 0);
    //fit.setlimit(3, 1.0e3, 1);
    //std::cout << "begin fitting" << std::endl;
    //fit.lmfit();
    fit.siman();
    //fit.info();
    //std::cout << fit.iternum << std::endl;
    //for(int i = 0; i < fit.npout; ++i) std::cout << fit.pout[i] << " ";
    //std::cout << std::endl;

    //line_gaussian line(4861);
    //line_lorentzian line2(5007);
    
    //func = &line;
    //func->info();
    //func->calc(num, x, y, p);
    //func->info();
    
    //func = &line2;
    //func->info();
    //func->calc(num, x, y, p);
    //func->info();
    //std::cout << func->npar << std::endl;
    //std::cout << func->profile << std::endl;
    
    std::exit(0);
    // plotting
    for (int i = 0; i < num; ++i)
    {
        plotx[i] = x1[i];
        ploty[i] = y1[i];
    }
    
    //std::exit(0);
    // cpgbeg(0, "?", 1, 1);
    // cpgsvp(0.15, 0.95, 0.15, 0.95);
    // cpgswin(1000.0, 8000.0, -0.02, 0.1);
    // cpgline(num, plotx, ploty);
    /*
    for (int i = 0; i < num; ++i)
    {
        plotx[i] = x2[i];
        ploty[i] = y2[i];
    }
    cpgsci(4);
    cpgline(num, plotx, ploty);*/
    // cpgsci(1);
    // cpgbox("bcnst", 0.0, 0, "bcnstv", 0.0, 0);
    // cpglab("x", "y", "");
    // cpgend();
    return 0;
}
