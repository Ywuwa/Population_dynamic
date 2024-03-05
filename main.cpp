#include <iostream>

#include "inout.h"
#include "clustering.h"
#include "mfunctions.h"
#include "selective_grad.h"
#include <random>
int main() {
    //================================== MODEL INITIALIZING ===========================================
    Model model;
    const std::string filename    = "model_param.txt";
    std::string some_prefix = "";
    std::string filename2   = some_prefix + ".txt";
    int err_code = 0;
    const auto start{std::chrono::steady_clock::now()}; // start point
    int num = omp_get_num_procs();
    printf("number of available cores: %d\n", num);
    printf("init strated\n");
    model.dim = 2;

    err_code = model_input(model, filename);
    if (err_code < 0) {
        printf("MODEL INPUT ERROR CODE: %d\n", err_code);
        printf("tap any number to continue: ");
        scanf("%d", &err_code);
        return 0;
    }

    unsigned int n = 1;
    const double b = 3.0;
    const double dx = 0.1;
    for (int i=0; i<model.dim; i++) n *= 60;

    unsigned int nstart = 0;
    unsigned int nend = 0;
    printf("input start file number and end file number (>0, //500): \n");
    scanf("%u %u", &nstart, &nend);
    std::string postfix = std::to_string(nstart);
    some_prefix = "res_dir/" + postfix;
    err_code = dense_input(model, some_prefix + filename2);
    if (err_code < 0) {
        printf("RESOURCE DENSITY INPUT ERROR CODE: %d\n", err_code);
        if (err_code==-2){
            printf("default input: R(y) = K(y)\n");
            std::unique_ptr<double []> a = std::make_unique<double []> (model.dim);
            model.R = std::make_unique<double []>(n);
            for (unsigned int j=0; j<n; j++){
                        a[0] = -b + (j/60)*dx;
                        a[1] = -b + (j%60)*dx;
                        model.R[j] = Kore(a.get(), 1.0, model.k0, model.dim);
            }
        }
        else{
            printf("tap any number to continue: ");
            scanf("%d", &err_code);
            return 0;
        }
    }
    printf("MODEL INITIALIZING IS FINISHED:\n");
    printf("//==  dimensions quantity: %d\n", model.dim+1);
    if (model.rtype) printf("//==  resource type: logistic resource\n");
    else printf("//==  resource type: chemostat resource\n");
    printf("//==  max carrying capacity: %d\n", model.k0);
    printf("//==  sig_a: %f, sig_g: %f\n", model.sig_a, model.sig_g);
    printf("//==  eps: %f, xi: %f\n", model.eps, model.xi);
    printf("//==  recreation rate: %f, dissipation rate: %f, death rate: %f\n", model.birth, model.dissip, model.death);
    printf("//==  m1: %f, m2: %f\n", model.m[0], model.m[1]);

    printf("//==  R[0][0] = %10.5e\n", model.R[0]);
    //=================================================================================================
    //================================== POPULATIONZ INITIALIZING =====================================
    std::vector<population> poplist;
    some_prefix = "populationz/" + postfix;
    std::cout << some_prefix << std::endl;
    err_code = populationz_input(poplist, model.dim, some_prefix + filename2);
    if (err_code < 0) {
        printf("POPULATIONZ INPUT ERROR CODE: %d\n", err_code);
        printf("tap any number to continue: ");
        scanf("%d", &err_code);
        return 0;
    }
    printf("POPULATIONZ INITIALIZING IS FINISHED:\n");
    printf("//==  total number of populations: %llu\n", poplist.size());
    for (size_t i = 0; i < poplist.size(); i++) {
        poplist[i].descrip();
    }
    info(poplist);
    //=================================================================================================
    printf("if u want to compute population dynamic print 0\n");
    printf("if u had already compute it and need to the selective gradient dynamic print 1\n");
    unsigned int var = 0;
    scanf("%u", &var);
    if (var==0)
    { printf("Population dynamic computing started\n");
    //======================================= GENERAL PART ============================================
    std::vector<population> buf_poplist;
    for (unsigned int i=nstart+1; i<nend+1; i++)
    {
        if (poplist.size()<1) {printf("ALL POPULATIONZ ARE EXTINCTED\n"); break;}
        //------------------------ calculate resource ------------------------------
        omp_set_num_threads(num);
        #pragma omp parallel
        {
            std::unique_ptr<double []> y = std::make_unique<double []> (model. dim);
            #pragma omp for schedule(dynamic,15)
            for (unsigned int j=0; j<n; j++){
                y[0] = -b + (j/60)*dx;
                y[1] = -b + (j%60)*dx;
                model.R[j] = runge_compute_res(y.get(), poplist, model.R[j], model, poplist.size());
            }
        }
        //-------------------------------------------------------------------------
        //---------------------- calculate populationz ----------------------------
        //buf_poplist = std::move(poplist);
        std::swap(buf_poplist, poplist);

        paral_runge_pop(buf_poplist, model.R.get(), model.m.get(),
                        model.xi, model.death, model.eps, model.sig_a, model.sig_g,
                        model.dim, buf_poplist.size());

        const size_t len = buf_poplist.size();
        for (size_t j=len; j>0; j--){
            if (buf_poplist[j-1].N > 1e-06) { poplist.push_back(std::move(buf_poplist[j-1])); buf_poplist.pop_back(); }
            else buf_poplist.pop_back();
        }
        //-------------------------------------------------------------------------
        if (i%500 == 0){
            pop_clust2(poplist, 0.014142);
            some_prefix = "res_dir/" + std::to_string(i);
            err_code = dense_output(model, some_prefix + filename2);
            if (err_code < 0) {
                printf("RESOURCE DENSITY OUTPUT ERROR CODE: %d\n", err_code);
                printf("ITERATION %d\n", i);
                printf("tap any number to continue: ");
                scanf("%d", &err_code);
                return 0;
            }
            some_prefix = "populationz/" + std::to_string(i);
            err_code = populationz_output(poplist, model.dim, some_prefix + filename2);
            if (err_code < 0) {
                printf("POPULATIONZ OUTPUT ERROR CODE: %d\n", err_code);
                printf("tap any number to continue: ");
                scanf("%d", &err_code);
                return 0;
            }
            printf("%u/%u iteration passed, current populationz amount: %llu\n", i, nend, poplist.size());
            info(poplist);
        }
        if (poplist.size()<700){
            const int choice = pop_choice(poplist);
            population p;
            p.mutate(poplist[choice]);
            poplist.push_back(std::move(p));
        }
    }
    //=================================================================================================
    }
    if (var==1)
    { printf("Selective gradient computing started\n");
    //=================================== SELECTIVE GRADIENT DYNAMIC ==================================
    std::vector<double> sgrad;
    for (unsigned int i=nstart; i<nend; i+=500)
    {
        poplist.clear();
        postfix = std::to_string(i);
        some_prefix = "res_dir/" + postfix;
        err_code = dense_input(model, some_prefix + filename2);
        if (err_code < 0) {
            printf("RESOURCE DENSITY INPUT ERROR CODE: %d\n", err_code);
            if (err_code==-2){
                printf("default input: R(y) = K(y)\n");
                std::unique_ptr<double []> a = std::make_unique<double []> (model.dim);
                model.R = std::make_unique<double []>(n);
                for (unsigned int j=0; j<n; j++){
                            a[0] = -b + (j/60)*dx;
                            a[1] = -b + (j%60)*dx;
                            model.R[j] = Kore(a.get(), 1.0, model.k0, model.dim);
                }
            }
            else{
                printf("tap any number to continue: ");
                scanf("%d", &err_code);
                return 0;
            }
        }
        some_prefix = "populationz/" + postfix;
        err_code = populationz_input(poplist, model.dim, some_prefix + filename2);
        if (err_code < 0) {
            printf("POPULATIONZ INPUT ERROR CODE: %d\n", err_code);
            printf("tap any number to continue: ");
            scanf("%d", &err_code);
            return 0;
        }
        const double d = selective_grad(model, poplist);
        printf("%u/%u. Populations amount: %llu. Selective gradient metric: %lf\n", i, nend, poplist.size(), d);
        sgrad.push_back(d);
    }
    some_prefix = "populationz/selective_gradient_dynamic";
    err_code = grad_output(sgrad, some_prefix + filename2);
    if (err_code < 0) {
                printf("GRADIENT METRIC OUTPUT ERROR CODE: %d\n", err_code);
                printf("tap any number to continue: ");
                scanf("%d", &err_code);
                return 0;
    }
    //=================================================================================================
    }
    //================================== ENFREE MEMORY & FINAL OUTPUT =================================

    printf("Output is completed\n");
    printf("Done\n");

    const auto end{std::chrono::steady_clock::now()}; // end point
    const std::chrono::duration<double> elapsed_seconds{end-start};

    printf("Total time: %lf. Input any number to continue: ", elapsed_seconds);
    scanf("%d", &err_code);
    return 0;
}
