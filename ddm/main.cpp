//
//  main.cpp
//  ddm
//
//  Created by Zhao Yongfeng on 14/11/3.
//  Copyright (c) 2014年 ZYF. All rights reserved.
//

#include "functions.h"

//The name prefix of the images should be given in bash as argv. (e.g.: IM_340_01_X)
//Or, give "simulation" can deal with the simulation data.
//Or, give "recover" can do the fit by reading datag.txt and q.txt directly, save doing time-consuming FFT and averaging.
int main(int argc, const char * argv[])
{
    omp_set_num_threads(8);     //Number of threads
    stringstream arg;		//To read the argv
    arg << argv[1];		 
    
//    arg << "simulation";
    gsl_matrix* datag;		//g(q,t) matrix. Need to claim at the beginning since we may need to deal with simulation data.
    vector<double> qabs;	//Absolute value of q array.
    qabs.reserve((dimk - 1)*(dimk - 1) / 2);	//the number of effective q will not exceed (dimk - 1)*(dimk - 1)/2.
    int qsize;				//Element number of q array.
    
    if (arg.str() == "recover")		//recover from pre-calculated g(q,t) data
    {
        cout << "Loading datag.txt..." << endl;
        
        ifstream qfile("q.txt");
        while (!qfile.eof())	//Reading q array
        {
            double tempq;
            qfile >> tempq;
            qabs.push_back(tempq);
        }
        qfile.close();
        qsize = (int)qabs.size();
        datag = gsl_matrix_alloc(qsize, numOfDiff);		//Allocate memory for g(q,t) matrix.
        
        ifstream datagfile("datag.txt");
        
        for (int iterq = 0; iterq < qsize; ++iterq)		//Reading g(q,t) matrix.
        {
            for (int itertau = 0; itertau < numOfDiff; ++itertau)
            {
                double tempdata;
                datagfile >> tempdata;
                gsl_matrix_set(datag, iterq, itertau, tempdata);
            }
        }
        datagfile.close();
    }
    else
    {
        string filePrefix = arg.str();
        //string filePrefix="cl-m1-ddm-01_X";
        //		vector<gsl_matrix*> imageSeq(numOfSeq);
        
        gsl_vector* aveVec = gsl_vector_alloc(dim);
        gsl_vector_set_all(aveVec, 1.0 / dim);			//Vecter used in calculating average. I hope BLAS can help speed up the average.
        
        //		cout << "Loading files..." << endl;
        
        
        //cout << "Shifting images..." << endl;
        //for (int iter = 0; iter < numOfDiff; ++iter)
        //{
        //	double max = 0;
        //	int drow = 0;
        //	int dcol = 0;
        //	for (int iterrow = -winDim; iterrow <= winDim; ++iterrow)
        //	{
        //		for (int itercol = -winDim; itercol <= winDim; ++itercol)
        //		{
        //			double corr = shiftCorrelation2d(imageSeq[iter], imageSeq[iter + 1], iterrow, itercol, aveVec);
        //			if (corr > max)
        //			{
        //				max = corr;
        //				drow = iterrow;
        //				dcol = itercol;
        //			}
        //		}
        //	}
        //	shiftImage(&imageSeq[iter + 1], drow, dcol);
        
        //	//        cout << "Shifting images..." << (iter+1.0)/numOfDiff*100 << "% finished" << endl;
        //}
        
        cout << "Loading and calculating FFT..." << endl;
        
        vector<gsl_matrix_complex*> imageSeqk(numOfSeq);	//Sequence for storing image after FFT.
        
#pragma omp parallel
        {
            gsl_matrix* fftMatrix = gsl_matrix_alloc(dim, dim);
            gsl_matrix_complex* resultMatrix = gsl_matrix_complex_alloc(dim, dimk);		
            //dimk=dim/2+1, FFTW only return around half points for r2c FFT to save memory usage.
            
            fftw_plan fft2plan;		//FFT need to make a plan before reading data.
#pragma omp critical (make_plan)
            fft2plan = fftw_plan_dft_r2c_2d(dim, dim, fftMatrix->data, (fftw_complex *)resultMatrix->data, FFTW_MEASURE);
            
#pragma omp for
            for (int iter = 0; iter < numOfSeq; ++iter)
            {
                imageSeqk[iter] = gsl_matrix_complex_calloc(dim, dimk);
                gsl_matrix* tempMatrix;     //Use temperary matrix to free the memory allocated in functions correctly.
                                            //Copy operation is not the bottleneck of the calculation, so it's not necessary to optimize this part.
                                            //But I/O is the bottleneck.
                if (filePrefix=="simulation")
                {
                    stringstream fileName;
                    fileName << iter << ".txt";
                    tempMatrix = readSim(fileName.str());
                }
                else
                {
                    stringstream fileName;
                    fileName << filePrefix << iter+1 << ".tif";		//Making path
                    tempMatrix = readTiff(fileName.str());
                }
                
                gsl_matrix_memcpy(fftMatrix, tempMatrix);
                gsl_matrix_free(tempMatrix);	//The problem is the allocation of matrix is inside the readSim/readTif function.
                
                fftw_execute(fft2plan);
                
                gsl_matrix_complex_memcpy(imageSeqk[iter], resultMatrix);
            }
            gsl_matrix_complex_free(resultMatrix);
            fftw_destroy_plan(fft2plan);
        }
        //		imageSeq.clear();
        
        cout << "Calculating average of square module for different tau... 0% finished." << endl;
        vector<gsl_matrix*> imagekDiff(numOfDiff);		//For storing the time difference of the imageSeqk
        int progress = 0;       //Indicator of the progess.
#pragma omp parallel for
        for (int iterdiff = 1; iterdiff <= numOfDiff; ++iterdiff)
        {
            gsl_matrix* temp = gsl_matrix_alloc(dim, dimk);
            gsl_matrix_set_zero(temp);
            for (int itert = 0; itert < numOfSeq - numOfDiff; ++itert)
            {
                for (int itermem = 0; itermem < numOfk; ++itermem)
                {
                    double* later = imageSeqk[itert + iterdiff]->data;		//I(q, t+\tau)
                    double* before = imageSeqk[itert]->data;				//I(q, t)
                    double real = later[itermem * 2] - before[itermem * 2];	//Difference of Re
                    double image = later[itermem * 2 + 1] - before[itermem * 2 + 1];	//Difference of Im
                    temp->data[itermem] += real*real + image*image;			//Noticed |a|=|a*|
                }
            }
            progress += 1; // numOfSeq - iterdiff;
            gsl_matrix_scale(temp, 1.0 / (numOfSeq - numOfDiff));		//Average on t
            imagekDiff[iterdiff - 1] = temp;
            cout << "Calculating average of square module for different tau... " << 100.0*progress / numOfDiff << "% finished." << endl;
        }
        
#pragma omp parallel for
        for (int iter = 0; iter < numOfSeq; ++iter)		//Free the memory
            gsl_matrix_complex_free(imageSeqk[iter]);
        imageSeqk.clear();
        
        cout << "Averaging on directions of q..." << endl;
//        vector<int> q;
//        q.reserve((dimk - 1)*(dimk - 1) / 2);
//        q.push_back(0);
//        for (int iteri = 1; iteri < dimk; ++iteri)
//        {
//            for (int iterj = 0; iterj <= iteri; ++iterj) {
//                int index = iteri*iteri + iterj*iterj;
//                if (find(q, index) == -1)
//                {
//                    q.push_back(index);
//                }
//            }
//        }
//        
//        qsize = (int)q.size();
//        //sorting q
//        int stackIter = 0;
//        int* stack = new int[qsize];
//        stack[stackIter] = 0;
//        ++stackIter;
//        stack[stackIter] = qsize - 1;
//        
//        while (stackIter > -1)
//        {
//            int right = stack[stackIter];
//            --stackIter;
//            int left = stack[stackIter];
//            --stackIter;
//            int empty = left - 1;
//            double pivotvalue = q[right];
//            for (int64_t iter = left; iter < right; ++iter)
//            {
//                if (q[iter] < pivotvalue)
//                {
//                    ++empty;
//                    int temp = q[iter];
//                    q[iter] = q[empty];
//                    q[empty] = temp;
//                }
//            }
//            int temp = q[right];
//            q[right] = q[empty + 1];
//            q[empty + 1] = temp;
//            if (left < empty)
//            {
//                ++stackIter;
//                stack[stackIter] = left;
//                ++stackIter;
//                stack[stackIter] = empty;
//            }
//            if (empty + 2 < right)
//            {
//                ++stackIter;
//                stack[stackIter] = empty + 2;
//                ++stackIter;
//                stack[stackIter] = right;
//            }
//        }
//        delete[] stack;
//        
//        datag = gsl_matrix_calloc(qsize, numOfDiff);
//        gsl_matrix* count = gsl_matrix_alloc(qsize, numOfDiff);
//        gsl_matrix_set_zero(datag);
//        gsl_matrix_set_zero(count);
//#pragma omp parallel for
//        for (int itertau = 0; itertau < numOfDiff; ++itertau)
//        {
//            for (int iterrow = 0; iterrow < dim; ++iterrow)
//            {
//                int kx;
//                if (iterrow < dimk)
//                    kx = iterrow;
//                else
//                    kx = dim - iterrow;
//                for (int itercol = 0; itercol < dimk; ++itercol)
//                {
//                    int ky = itercol;
//                    int factor = 2;
//                    if (ky % (dimk - 1) == 0)
//                        factor = 1;
//                    int k2 = kx*kx + ky*ky;
//                    int index = quickFind(q, k2);
//                    
//                    double temp = gsl_matrix_get(imagekDiff[itertau], iterrow, itercol)*factor;
//                    temp += gsl_matrix_get(datag, index, itertau);
//                    gsl_matrix_set(datag, index, itertau, temp);
//                    double tempindex = gsl_matrix_get(count, index, itertau) + factor;
//                    gsl_matrix_set(count, index, itertau, tempindex);
//                }
//            }
//        }
        
        double qmax=sqrt(dim*dim/2)*dq;
        qsize=floor(qmax/qstep)+1;
        qabs.resize(qsize);     //qabs is the absolute value of q
        for (int iter=0; iter<qsize; ++iter)
        {
            qabs[iter]=(iter+0.5)*qstep;    //The midpoint of the cirque
        }
        
        datag = gsl_matrix_calloc(qsize, numOfDiff);
        gsl_matrix* count = gsl_matrix_alloc(qsize, numOfDiff);
        gsl_matrix_set_zero(datag);
        gsl_matrix_set_zero(count);     //Number of elements
#pragma omp parallel for
        for (int itertau = 0; itertau < numOfDiff; ++itertau)
        {
            for (int iterrow = 0; iterrow < dim; ++iterrow)
            {
                int kx;
                if (iterrow < dimk)		//kx=[0, 1, ..., dim/2, 1-dim/2, 2-dim/2, ..., -1]
                    kx = iterrow;
                else
                    kx = iterrow - dim;
                for (int itercol = 0; itercol < dimk; ++itercol)
                {
                    int ky = itercol;
                    int factor = 2;     //If the data is not the q=0 or N-1, it appears twice in q-space for a_i=a_{dim-i}.
                    if (ky % (dimk - 1) == 0)
                        factor = 1;
                    double k2 = kx*kx + ky*ky;
                    int index = floor(sqrt(k2)*dq/qstep);
                    
                    //This is fast enough. No need to speed up.
                    double temp = gsl_matrix_get(imagekDiff[itertau], iterrow, itercol)*factor;
                    temp += gsl_matrix_get(datag, index, itertau);
                    gsl_matrix_set(datag, index, itertau, temp);		//datag[index, itertau]+=imagekDiff[itertau][tierrow, itercol]*factor
                    double tempcount = gsl_matrix_get(count, index, itertau) + factor;
                    gsl_matrix_set(count, index, itertau, tempcount);	//count[index, itertau]+=factor
                }
            }
        }
        gsl_matrix_div_elements(datag, count);		//Average
#pragma omp parallel for
        for (int iter = 0; iter < numOfDiff; ++iter)
            gsl_matrix_free(imagekDiff[iter]);
        imagekDiff.clear();
        
        cout << "Printing unfitted data..." << endl;
        ofstream qfile("q.txt");
        for (int iter = 0; iter < qsize; ++iter)
        {
//            qfile << sqrt(q[iter])*dq << endl;
            qfile << qabs[iter] << endl;
        }
        qfile.close();
        ofstream datagfile("datag.txt");
        for (int iterq = 0; iterq < qsize; ++iterq)
        {
            for (int itertau = 0; itertau < numOfDiff; ++itertau)
            {
                datagfile << gsl_matrix_get(datag, iterq, itertau) << " ";
            }
            datagfile << endl;
        }
        datagfile.close();
        //        qabs.resize(qsize);
        //        for (int iterq = 0; iterq < qsize; ++iterq)
        //        {
        ////            qabs[iterq] = sqrt(q[iterq])*dq;
        //        }
        gsl_vector_free(aveVec);
    }
    
#ifdef NeedLaplaceTrans
    cout << "Numerical Laplace transforming..." << endl;
    gsl_matrix* transM=gsl_matrix_alloc(numOfDiff, numOfDiff);		//Transformation matrix
    for (int itert = 0; itert < numOfDiff; ++itert)
    {
        for (int iters = 0; iters < numOfDiff; ++iters)
        {
            double tau=(itert+1)*dt;
            double s=exp(-1+0.001*iters);		//s is sampled in log scale
            gsl_matrix_set(transM, itert, iters, 0.01*exp(-tau*s));		//\Delta t=0.01, numerical integration here.
        }
    }

    ofstream sfile("s.txt");
    for (int iters = 0; iters < numOfDiff; ++iters)
    {
        double s=exp(-1+0.001*iters);
        sfile << s << endl;
        gsl_matrix_set(transM, numOfDiff-1, iters, 0.5*gsl_matrix_get(transM, numOfDiff-1, iters));
        gsl_matrix_set(transM, 0, iters, 1.5*gsl_matrix_get(transM, 0, iters));		//0.5+1. 1 is the estimation of the contribution in tau=(0,0.01], which is not computable since tau>0.
    }
    sfile.close();

    gsl_matrix* temp=gsl_matrix_alloc(qsize, numOfDiff);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, datag, transM, 0, temp);		//LT is performed by matrix multiplication. Call BLAS to calculate matrix multiplication.
    gsl_matrix_free(transM);

    ofstream datagfile("ldatag.txt");		//Print Laplace transformed g(q,s)
    for (int iterq = 0; iterq < qsize; ++iterq)
    {
        for (int itertau = 0; itertau < numOfDiff; ++itertau)
        {
            datagfile << gsl_matrix_get(temp, iterq, itertau) << " ";
        }
        datagfile << endl;
    }
    datagfile.close();

#endif

    cout << "Fitting..." << endl;
    
    const gsl_multifit_fdfsolver_type *solverType;	//GSL solver
    solverType = gsl_multifit_fdfsolver_lmsder;		
    //Using Levenberg-Marquardt algorithm as implemented in the scaled lmder routine in minpack. Jacobian is given.
    
    gsl_matrix* fittedPara=gsl_matrix_alloc(qsize-1, numOfPara);	//To store the fitting result and error.
    gsl_matrix* fitErr=gsl_matrix_alloc(qsize-1, numOfPara);
    int* status = new int[qsize-1];		//Record the status of fitting.
    
    double tau[num_fit];
    for (int itertau=0; itertau<num_fit; ++itertau)
    {
        tau[itertau]=(itertau+1)*dt;
#ifdef ISFRUNANDTUMBLE_3D
        tau[itertau]=exp(-1+0.001*itertau);
#endif
#ifdef ISFRunAndTumbleAndDiffusion
        tau[itertau]=exp(-1+0.001*itertau);
#endif
    }

    //Initial guess
#ifdef ISFSWIMMER
    double inipara[numOfPara]={0.7, 0.4, 5, 2, 2e12, 1e10};
#endif
#ifdef ISFSWIMMERSIMPLER
    double inipara[numOfPara]={0.5, 0.4, 10, 1, 1};
#endif
#ifdef ISFRUNANDTUMBLE
    double inipara[numOfPara] = {1, 10, 2e12, 1e10 };
#endif
#ifdef ISFRUNANDTUMBLE_3D
    double inipara[numOfPara] = {5, 0.5, 2e12, 1e10 };
#endif
#ifdef ISFRunAndTumbleAndDiffusion
    double inipara[numOfPara] = {0.5, 5, 0.5, 1, 2e12, 1e10 };
#endif
    
    //Get the selected data for fitting
    gsl_matrix* datafit = gsl_matrix_alloc(qsize, num_fit);
    for (int iterq = 0; iterq < qsize; ++iterq)
    {
        for (int iterf = 0; iterf < num_fit; ++iterf)
        {
            gsl_matrix_set(datafit, iterq, iterf, log(gsl_matrix_get(temp, iterq, iterf)));		//Fitting in log scale.
        }
    }
    int progress=0;		//Indicator of progress.
    //	ofstream debugfile("debug.txt");
#pragma omp parallel for
    for (int iterq=1; iterq<qsize; ++iterq)
    {
        gsl_multifit_function_fdf fitfun;		//Function point.
        gsl_vector_view dataAry=gsl_matrix_row(datafit, iterq);
        dataStruct sdata;		//GSL data structure
        sdata.data=dataAry.vector.data;
        sdata.tau=tau;
        sdata.q=qabs[iterq];
        
        //API
        fitfun.f=&ISFfun;
        fitfun.df=&dISFfun;
        fitfun.fdf=&fdISFfun;
        fitfun.n=num_fit;
        fitfun.p=numOfPara;
        fitfun.params=&sdata;
        
        //Estimation of A(q) and B(q)
	    inipara[5] = gsl_matrix_get(datag, iterq, 0);
        inipara[4] = gsl_matrix_get(datag, iterq, 3000)-inipara[5];
        
        //Initiallization of the solver
        gsl_vector_view para=gsl_vector_view_array(inipara, numOfPara);
        gsl_multifit_fdfsolver* solver = gsl_multifit_fdfsolver_alloc(solverType, num_fit, numOfPara);
        gsl_multifit_fdfsolver_set(solver, &fitfun, &para.vector);
        int iter=0;
        //gsl_vector* g=gsl_vector_alloc(numOfPara);
        do
        {
            gsl_multifit_fdfsolver_iterate(solver);		//Iterate one step.
            status[iterq-1] = norm0_rel_test(solver->dx, solver->x, 1e-10, 1e-10);		//Test the exiting condition

            //gsl_multifit_gradient(solver->J,solver->f, g);
            //status[iterq-1]=gsl_multifit_test_gradient(g, 1e-5);
//			status[iterq - 1] = covar_rel_test(solver->J, solver->x, 1e-4);
            
            ++iter;
            if (iter>maxIter)
            {
                status[iterq-1]=GSL_EMAXITER;
            }
        } while (status[iterq-1] == GSL_CONTINUE);
        //gsl_vector_free(g);

        //Estimating the error.
        gsl_matrix* covar=gsl_matrix_alloc(numOfPara, numOfPara);
        gsl_multifit_covar(solver->J, 0.0, covar);
        for (int iterpara=0; iterpara<numOfPara; ++iterpara)	//Record result.
        {
            gsl_matrix_set(fittedPara, iterq-1, iterpara, gsl_vector_get(solver->x, iterpara) );
            gsl_matrix_set(fitErr, iterq-1, iterpara, sqrt(gsl_matrix_get(covar, iterpara, iterpara)) );    //Not presice in log scale
        }
        gsl_matrix_free(covar);
        gsl_multifit_fdfsolver_free(solver);
        
        progress+=1;
        cout << "Fitted q=" << qabs[iterq] << " at iter=" << iter << ", " << 100.0*progress / qsize << "% completed from core No." << omp_get_thread_num() << ", "<< gsl_strerror(status[iterq-1]) << "." << endl;
    }
    
    //    cout << "Printing fit result..." << endl;
    ofstream fitparafile("fitparafile.txt");
    ofstream fiterrfile("fiterrfile.txt");
    ofstream statusfile("status.txt");
    for (int iterq=0; iterq<qsize-1; ++iterq)
    {
        for (int iterpara=0; iterpara<numOfPara; ++iterpara)
        {
            fitparafile << gsl_matrix_get(fittedPara, iterq, iterpara) << " ";
            fiterrfile << gsl_matrix_get(fitErr, iterq, iterpara) << " ";
        }
        fitparafile << endl;
        fiterrfile << endl;
        statusfile << iterq << ": q=" << qabs[iterq+1] << ", "<< gsl_strerror(status[iterq]) << endl;
    }
    
    fitparafile.close();
    fiterrfile.close();
    
    delete[] status;
    gsl_matrix_free(temp);
    gsl_matrix_free(datag);
    gsl_matrix_free(fittedPara);
    gsl_matrix_free(fitErr);
    gsl_matrix_free(datafit);
    return 0;
}