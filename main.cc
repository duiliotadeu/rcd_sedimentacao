/* 
 * OLD VERDION 
#include "config.h"

#include "RCDUtil.h"
#include "RiemannProblem.h"
#include "bdry.h"
#include "dct.h"
#include "keymatch.h"
#include "newton.h"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
*/

// 
// nova versao do modelo de sedimentacao
//

// #include "config.h"

#include "RCDUtil.h"
#include "RiemannProblem.h"
#include "bdry.h"
#include "dct.h"
#include "keymatch.h"
#include "newton.h"
#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <limits>

void Verifica(std::string foldername)
{
    // Subrotina para verificar se o subdiretorio para gravação existe
    std::string filename(foldername);
    filename.append("/check.txt");
    std::cout << "Filename = " << filename << std::endl;
    std::ofstream myfile;
    myfile.open(filename);
    if (myfile.fail())
    {
        std::cerr << "[Abortando]  Subdiretorio " << foldername << " nao existe" << std::endl;
        std::cerr << "Crie o subdiretorio para poder executar." << std::endl;
        exit(1);
    }
    else
    {
        myfile.close();
        const char * mystr = filename.c_str();
        remove(mystr); // remove arquivo que foi criado para testar.
        std::cout << "Subdiretorio existe" << std::endl;
    }
}



class Sedimentacao : public PhysicsRCD
{
public:
    void setConstantValues(PhysicsDataHolder& phyDH) override;
    
    void jet(PhysicsDataHolder& phyDH, double* Vj, double* Un, double& p, double& u,
            const double& x, const double& t, const std::size_t& idx) override;
            
    void setPDEAndCstFunctions(PhysicsDataHolder& phyDH, double* Un, double& p,
            double& u, const double& x, const double& t,
            const std::size_t& idx) override;

    /*
     * A seguir metodos e dados especificos do problema. O  RCD não exige estas
     *declaracoes,
     * elas facilitam a implementacao das funcoes acima.
     */

    /*
     * Constantes do problema
     */
    double *diam;  // particle diameter
    double *rho_p; // particle density
    double rho_f;  // Fluid density
    double mu_f;   // fluid viscodity
    double settling_n; // Hindered Settling Function : parameter n
    double settling_lambda; // Hindered Settling Function : paramenter lambda
    int Ne;  // Number of particle species
    // Parametros para definir a função de fluxo
    double M;   //  = - g d_1^2/(18*mu_f)
    double *delta; //  delta_i = d_i^2 / d_1^2
    double tol_numdiffusion;  // Numerical Diffusion 
    double phi_max; // maximum total concentration

    /*
     * Construtor da classe
     */
    // Primeiro, os parametros que pertence a qualquer fisica trifasica
    Sedimentacao(
            double *diam_, double *rho_p_, double rho_f_, double mu_f_,
            double settling_n_, double settling_lambda_, int Ne_, 
            double tol_numdiffusion_, double phi_max_)
        : PhysicsRCD()
        ,
        rho_f(rho_f_), mu_f(mu_f_), settling_n(settling_n_),
        settling_lambda(settling_lambda_), Ne(Ne_), tol_numdiffusion(tol_numdiffusion_),
        phi_max(phi_max_)
    {
        // ALOCACAO DE DIAM  E  RHO       
        diam = new double[Ne_];
        rho_p = new double[Ne_];
        delta = new double[Ne_];
    
        M = - 9.81*diam_[0]*diam_[0]/(18*mu_f_);
        for (int i=0; i<Ne_; i++)
        {
            diam[i] = diam_[i];
            rho_p[i] = rho_p_[i];
            delta[i] = diam_[i]/diam_[0];
        }  
    }

    ~Sedimentacao()
    {
        delete [] diam;
        delete [] rho_p;
        delete [] delta;
    }
};


void Sedimentacao::setConstantValues(PhysicsDataHolder& phyDH)
{
	// Introducao de difusao numerica no modelo
	// O valor é definido no arquivo de entrada
	int i;
	
	for (i=0; i< (this->Ne * this->Ne) ;i++ )
	{
		phyDH.B[i]=0.0;
	}
	for (i=0; i<this->Ne ; i++)
		phyDH.B[i*this->Ne + i] = this->tol_numdiffusion;
	
    
    

	
    // Jacobiano da funcao de acumulacao, eh a identidade, nao colocamos os
    // valores nulos.
    for (i = 0; i<this->Ne; i++)
		phyDH.DG[i*this->Ne +i ] = 1.0;

}

void
Sedimentacao::jet(PhysicsDataHolder& phyDH, double* Vj, double* /* Un */,
        double& /* p */, double& /* u */, const double& /* x */,
        const double& /* t */, const std::size_t& /* idx */)
{
    /* OBS: Deve-se sempre usar Vj dentro dessa funcao, a menos que se queira
       usar um box scheme,
     * neste caso converse com Ismael.
     *
     * Vj representa a solucao que esta sendo construida atraves do metodo de
     *Newton. O metodo toma
     * como chute inicial o vetor Un, e ao fim de sua execucao, aproxima o vetor
     *Un+1. A funcao jet
     * deve preencher todos os termos que dependem de Vj, isto eh,  os termos
     *que variam durante as
     * iteracoes do metodo de Newton.
     */
    
    double phi_tot {0.0}; // soma de todas as concentracoees
    double hsfunction {0.0};  // Valor da Hindered Settling Function
    double rho_mix {0.0};  // densidade da mistura
    int i, j;
    double A {0.0}; // = \sum_k delta[k]*phi[k]*(rho_p[k] - rho_mix)
    double tmp1, tmp2;
    double DV_Dphi;



    for (i=0;i<this->Ne;i++)
        phi_tot = phi_tot + Vj[i];

    rho_mix = (1-phi_tot)*this->rho_f;
    for (i=0;i<this->Ne;i++)
        rho_mix = rho_mix + this->rho_p[i]*Vj[i];

    // Evaluate Hindered Settling Function and Derivative
    if ((phi_tot < this->phi_max) && (phi_tot>0.0) )
    {
        hsfunction = pow( 1-phi_tot/this->phi_max  , this->settling_n);
        DV_Dphi = -(this->settling_n/this->phi_max) * pow(1-phi_tot/this->phi_max, this->settling_n - 1);
    }
    else
    {
        hsfunction = 0.0;
        DV_Dphi = 0.0; 
    }

    for (i=0; i<this->Ne; i++)
    {
        A = A + this->delta[i]*Vj[i]*(this->rho_p[i]-rho_mix);
    }

    tmp2 = 0.0;
    for (i=0; i<this->Ne; i++){
        tmp2 += this->delta[i] * Vj[i];
    }

    
    for (i=0; i<this->Ne; i++)
    {
        phyDH.F[i] = this->delta[i] * (this->rho_p[i]-rho_mix) - A; 
        phyDH.F[i] = this->M * hsfunction*Vj[i]*phyDH.F[i];
        
    }


    // Jacobiano do Fluxo
    for (i=0; i<this->Ne; i++)
    {
        tmp1 = this->delta[i]*(this->rho_p[i]-rho_mix) - A;
        
        for (j=0; j<this->Ne; j++)
        {
        // derivada D F_i / D phi_j
            phyDH.DF[i*this->Ne+j] = tmp1 * this->M * DV_Dphi * Vj[i] ; 
        }

        // Termo DF_i / D phi_i
        phyDH.DF[i*this->Ne+i] += tmp1 * this->M * hsfunction; 

        // Segunda parte da Derivada:
        for (j=0; j<this->Ne; j++) {
            phyDH.DF[i*this->Ne + j] += this->M * hsfunction * Vj[i] * 
                        ( -this->delta[i] * (this->rho_p[j] - this->rho_f) 
                          -this->delta[j] * (this->rho_p[j] - rho_mix)
                          + tmp2 * (this->rho_p[j]-this->rho_f)   ) ;
        }
    }
    
    // Termo G : Derivada do tempo
    for (i=0; i<this->Ne ; i++)
        phyDH.G[i] = Vj[i];

}

/*
 * A função setPDEAndCstFunctions é responsável por preencher os membros
 * de uma EDP, tais como acumulação, fluxo, difusão e reação
 */
void
Sedimentacao::setPDEAndCstFunctions(PhysicsDataHolder& phyDH, double* Un,
        double& /* p */, double& /* u */, const double& /* x */,
        const double& /* t */, const std::size_t& /* idx */)
{
    /*
     * Un representa a solucao que se tem como chute inicial no metodo de
     *Newton. Eh um vetor que
     * nao varia durante as iteracoes do metodo de Newton. Note que, diferente
     *da funcao jet, que
     * contem ambos Vj e Un,  esta funcao contem somente Un.   Isso porque essa
     *funcao eh somente
     * usada para preencher o vetor nao variante durante as iteracoes do metodo
     *de Newton.  O uso
     * correto dessa funcao e da funcao jet permite usar esquemas bastante
     *genrricos, como o box-
     * -scheme por exemplo (vide manual).
     */

    double phi_tot {0.0}; // soma de todas as concentracoees
    double hsfunction {0.0};  // Valor da Hindered Settling Function
    double rho_mix {0.0};  // densidade da mistura
    int i;
    double A {0.0}; // = \sum_k delta[k]*phi[k]*(rho_p[k] - rho_mix)


    for (i=0;i<this->Ne;i++)
        phi_tot = phi_tot + Un[i];

    rho_mix = (1-phi_tot)*this->rho_f;
    for (i=0;i<this->Ne;i++)
        rho_mix = rho_mix + this->rho_p[i]*Un[i];

    // Evaluate Hindered Settling Function and Derivative
    if ( (phi_tot < this->phi_max) && (phi_tot>0.0) )
    {
        hsfunction = pow( 1-phi_tot/this->phi_max  , this->settling_n);
        //DV_Dphi = -(n/this->phi_max) * pow(1-phitot/this->phi_max, this->settling_n - 1);
    }
    else
    {
        hsfunction = 0.0;
        //DV_Dphi = 0.0; 
    }

    for (i=0; i<this->Ne; i++)
    {
        A = A + this->delta[i]*Un[i]*(this->rho_p[i]-rho_mix);
    }

    for (i=0; i<this->Ne; i++)
    {
        phyDH.F[i] = this->delta[i] * (this->rho_p[i]-rho_mix) - A; 
        phyDH.F[i] = this->M * hsfunction*Un[i]*phyDH.F[i];
    }


    for (i=0; i<this->Ne; i++)
        phyDH.G[i] = Un[i];
}

/*
 * INITIAL SOLUTION
 */
void CalculateInitialSolution(double* V, std::size_t XSIZE, 
	int Ne, double *initial_concentration)
{
   size_t i,k;   
   for (i=0; i<XSIZE; i++)
   {
	   for (k=0;k<Ne;k++)
			V[i*Ne+k] = initial_concentration[k];
	}

}


void reorder_cotas(int num_cotapos, double *cota_position)
{
	double temp;
	for (int i=0; i<num_cotapos; i++)
		for (int j=i+1; j<num_cotapos; j++)
		{
			if (cota_position[i]>cota_position[j])
			{
				temp = cota_position[i];
				cota_position[i]=cota_position[j];
				cota_position[j]=temp;
			
			}
		}
}

void Set_idxcotas(double *x, std::size_t *idx_cotas, double *cota_position, std::size_t XSIZE, int num_cotapos)
{
	size_t pos{0}; 
	size_t j;
	int i;
	
	for (i=0; i<num_cotapos; i++)
	{
		idx_cotas[i]=0;
		for (j=pos; j<XSIZE; j++)
		{
			if (x[j] >= cota_position[i]) {
				std::cout << "Cota Idx::: " << j << " . " << x[j] << std::endl;
				idx_cotas[i]=j;
				pos = j+1;
				break;
			}
		} 
	}
}

/*
 * SAVE_SOLUTION_AT_TIME T
 * 
 * */
void Save_Solution(std::string dir, double *x, double *V, std::size_t XSIZE, int Ne, std::size_t timeidx, double current_time)
{
	std::string filename;
    std::string originaldir;
    originaldir = dir;
	// IMPORTANT: THIS VERSION IS FOR LINUX ONLY: NEED DIRECTIVE FOR WINDOWS
	filename = (dir.append("/SOLUTIONX_T")).append(std::to_string(timeidx));
	filename = filename.append(".out");
	std::cout << "SAVING AT TIME : " << current_time << " File : " << filename << std::endl;
	
    /* Impressao da solucao final:  */
    std::ofstream myfile(filename);
    myfile.precision(16);
    myfile << std::scientific;
    myfile << "Time: " << current_time << std::endl;
    myfile << "Nx: " << XSIZE << std::endl;
    myfile << "Species: " << Ne << std::endl;
    for (std::size_t i = 0, j = 0; i != XSIZE; i++, j += Ne) {
        myfile << x[i] << "   ";

        for (std::size_t k = 0; k != Ne; ++k)
            myfile << V[j + k] << "   ";
        myfile << std::endl;
    }
    myfile.close();

    
	filename = (originaldir.append("/DERSOLUTIONX_T")).append(std::to_string(timeidx));
	filename = filename.append(".out");

    std::cout << dir << std::endl;
    std::cout << "filename: >> "<< filename << std::endl;
    /* Impressao da Derivada solucao final:  */
    std::ofstream myfile2(filename);
    myfile2.precision(16);
    myfile2 << std::scientific;
    myfile2 << "Time: " << current_time << std::endl;
    myfile2 << "Nx: " << XSIZE << std::endl;
    myfile2 << "Species: " << Ne << std::endl;

    double sum1,sum2,dx;
    int idx;
    for (std::size_t i = 0; i != (XSIZE-1); i++) {
        myfile2 << x[i] << "   ";

        sum1=0.0;
        sum2=0.0;
        
        idx = i*Ne;
        for (std::size_t k=0;k<Ne; k++)
        {
            sum1 += V[idx+k];
            sum2 += V[idx+Ne+k];
        }
        myfile2 << (sum1-sum2)/dx << std::endl;

    }
    myfile2.close();
}

void scansave_topo(std::fstream &arquivo_topo, double currentTime, double *x, double *V, std::size_t XSIZE, int Ne)
{
    double sum1 = 0.0;
    double sum2 = 0.0;
    double tol_topo = 1.0;
    double numdiff1 = 0.0;
    double numdiff2 = 0.0;
    double numdiff3 = 0.0;

    double maxnumdiff=0.0;

    std::size_t position;
    std::size_t i,j,k, idx1, idx2;

    // Search for local maximum of modulus of num. derivative 
 /*   if (currentTime>0.0) 
    {
        idx1=(XSIZE-1)*Ne;
        idx2=(XSIZE-2)*Ne;
        for (k=0;k!=Ne;k++) {
            sum1 += V[idx1+k];
            sum2 += V[idx2+k];
        }
        numdiff1 = fabs(sum2-sum1);
        sum1=sum2;
        sum2=0.0;
        idx2=(XSIZE-3)*Ne;
        for (k=0;k!=Ne;k++) {
            sum2 += V[idx2+k];
        }
        numdiff2 = fabs(sum2-sum1);

        sum1 = sum2;
        sum2 = 0.0;
        for (i=(XSIZE-4); i>=0; i--)
        {            
            idx2=i*Ne;
            for (k=0;k!=Ne;k++) {
                sum2 += V[idx2+k];
            }
            numdiff3 = fabs(sum2 - sum1) ;
            sum1 = sum2;
            sum2 = 0.0;
            // arquivo_topo << currentTime << " " << x[i] << " " << numdiff1 << " " << numdiff2 << " " << numdiff3 << std::endl;
             if ((numdiff2 > numdiff1) && (numdiff2 > numdiff3)) 
            {
                 arquivo_topo << currentTime << " " << x[i]  << std::endl;
                 break;
            } 

            numdiff1=numdiff2;
            numdiff2=numdiff3;
            //std::cout << "x="<< x[i] << "  Num. Diff : " << numdiff1 << std::endl;
        }
        
    }
    else
    {
        double t=0.0;
        arquivo_topo << t << " " << x[XSIZE-1] << std::endl; 
    }
*/
    // NOVO: Media
    int idx_a, idx_b;
    double x_a, x_b;
    if (currentTime>0.0)
    {
        // calcula posicao acima da interface
        for (i=(XSIZE-1); i>=0; i--)
        {
            sum1=0.0;
            idx1=i*Ne;
            for (k=0; k<Ne; k++)
                sum1 +=  V[idx1+k];
            if (sum1 > 1.e-10)
            {
                idx_a = i;
                x_a = x[i];

                break;
            }
        }
        // calcula posicao abaixo da interface
        for (j=i; j>=0; j--)
        {
            sum1 = 0.0;
            idx1=j*Ne;
            for (k=0; k<Ne; k++)
                sum1 +=  V[idx1+k];
            sum2 = 0.0;
            idx2 = (j-1)*Ne;
            for (k=0;k<Ne;k++)
                sum2 += V[idx2+k];

            if (fabs(sum1-sum2) < 1.e-8)
            {
                idx_b = j;
                x_b = x[j];
                break;
            }
        }
        std::cout << "[TOPO] T="<<currentTime <<" x_a = " << x_a << "    x_b = "<< x_b << std::endl;
        arquivo_topo << currentTime << " " << ((x_a+x_b)*0.5) << std::endl;
        exit;
    }
    else
    {
        double t=0.0;
        arquivo_topo << t << " " << x[XSIZE-1] << std::endl;
    }

}


int main(int argc, char** argv)
{
   std::size_t NMAXITER;
   std::size_t PDEDim, CSTDim{0};
   std::size_t XSIZE, NOUTPUTS;
   //double MOVIETIME;
   double LENGTH;
   double DELTAT;
   std::size_t TSIZE;

   double XLEFT, XRIGHT, TFINAL;
   double ZEROTOL, DIFFTOL;
   
   double *initial_concentrations;
   double *diameters;
   double *rho_p;
   double rho_f;
   double mu_f;
   double settling_n, settling_lambda;
   int Ne;
   double tol_numdiffusion { 1.0e-5};
   double phi_max; 
   int num_cotapos;
   double *cota_position;
   std::size_t *idx_cotas;
   double *cota_time;
   int current_cota {0};
   int num_cotatimes;
   
   std::string foldername;
   
   /*
     * Class to read parameters from file
     */
   KeyMatch myKeyMatch;

   /*
     * From arguments provided from command line, decides which file use
     * as input file for the parameters below
     */
   myKeyMatch.define_in_out(argc, argv);

   myKeyMatch.define_in_out(argc, argv);
   myKeyMatch.register_entry(PDEDim, "*PDEDIM");
   // myKeyMatch.register_entry(CSTDim, "*CSTDIM");
   myKeyMatch.register_entry(XSIZE, "*XSIZE");
   myKeyMatch.register_entry(NOUTPUTS, "*NOUTPUTS"); // OK

//   myKeyMatch.register_entry(TSIZE, "*TSIZE"); // remover
//   myKeyMatch.register_entry(MOVIETIME, "*MOVIETIME");
   
   myKeyMatch.register_entry(LENGTH, "*LENGTH");
   myKeyMatch.register_entry(TFINAL, "*TFINAL");
   myKeyMatch.register_entry(DELTAT, "*DELTAT");
   
   myKeyMatch.register_entry(NMAXITER, "*NMAXITER");
   myKeyMatch.register_entry(ZEROTOL, "*ZEROTOL");
   myKeyMatch.register_entry(DIFFTOL, "*DIFFTOL");
   myKeyMatch.register_entry(foldername, "*OUTPUTFOLDER");
   myKeyMatch.register_entry(num_cotapos, "*NUMCOTAS");

   
   myKeyMatch.make();

   Ne = static_cast<int>(PDEDim);
   diameters = new double[Ne];
   rho_p = new double[Ne];
   initial_concentrations = new double[Ne];
   cota_position = new double[num_cotapos];
   idx_cotas = new std::size_t[num_cotapos];
   cota_time = new double[num_cotatimes];
   
   myKeyMatch.register_entry(diameters, PDEDim, "*DIAMETER");
   myKeyMatch.register_entry(rho_p, PDEDim, "*RHO_P");
   myKeyMatch.register_entry(initial_concentrations, PDEDim, "*PHI_INI");
   myKeyMatch.register_entry(mu_f, "*MU_FLUID");
   myKeyMatch.register_entry(rho_f, "*RHO_FLUID");
   myKeyMatch.register_entry(settling_n, "*HSFUNCTION_N");
   myKeyMatch.register_entry(settling_lambda, "*HSFUNCTION_LAMBDA");
   myKeyMatch.register_entry(phi_max, "*PHI_MAX");
   myKeyMatch.register_entry(tol_numdiffusion, "*TOL_NUMDIFF");
   myKeyMatch.register_entry(cota_position, num_cotapos, "*COTAS");

   myKeyMatch.make();


   //foldername="caso_a_8particulas";
// CHECK FOR FOLDER TO SAVE:
   std::cout << "Verifica folder"<< std::endl;
   Verifica(foldername);


   std::cout << "PDEDim    --> " << PDEDim << std::endl;
   std::cout << "CSTDim    --> " << CSTDim << std::endl;
   std::cout << "XSIZE(LIDO)     --> " << XSIZE << std::endl;
   std::cout << "NOUTPUTS  --> " << NOUTPUTS << std::endl;
   std::cout << "TSIZE     --> " << TSIZE << std::endl;
   std::cout << "XLEFT     --> " << XLEFT << std::endl;
   std::cout << "XRIGHT    --> " << XRIGHT << std::endl;
   std::cout << "TFINAL    --> " << TFINAL << std::endl;
   
   // Print Data
   std::cout << "Diameters  : " << std::endl;
   for (int i=0;i<Ne; i++)
      std::cout << "  Diam["<< i << "] = " << diameters[i] << std::endl;
   for (int i=0;i<Ne; i++)
      std::cout << "Particle Density ["<< i << "] = " << rho_p[i] << std::endl;
   for (int i=0;i<Ne; i++)
      std::cout << "Initial Concentration of ["<< i << "] = " << initial_concentrations[i] << std::endl;


   std::cout << "Fluid Viscosity : " << mu_f << std::endl;
   std::cout << "Fluid Density : " << rho_f << std::endl;
   std::cout << "H.S. Function - parameter n = "<< settling_n << std::endl;
   std::cout << "H.S.Function - Parameter Lambda = "<< settling_lambda << std::endl;
   std::cout << "Phi Maximum  = " << phi_max << std::endl;
   std::cout << "Numerical Diffusion Factor = " << tol_numdiffusion << std::endl;
   std::cout << " COTAS:  Quantity of X positions : " << num_cotapos << std::endl;
   reorder_cotas(num_cotapos, cota_position);
   for (int i=0; i< num_cotapos; i++)
		std::cout << "** Cota ["<< i<<"] = "<< cota_position[i] << std::endl;


   // definicao de XLEFT e XRIGHT
   XLEFT  = 0.0;
   XRIGHT = LENGTH;
   const double stpForMovie = NOUTPUTS > 0 ? TFINAL / (static_cast<int>(NOUTPUTS)) : -1.0;
   const double dx = ((XRIGHT - XLEFT) / (static_cast<double>(XSIZE)) ) > 0 ?
            ((XRIGHT - XLEFT) / (static_cast<double>(XSIZE)) ) :
            0;
   ++XSIZE;


   TSIZE = static_cast<size_t>( ceil( static_cast<double>(TFINAL)/static_cast<double>(DELTAT) ) );

   //tol_numdiffusion = 0.00002*dx;
   tol_numdiffusion  = 0.001*dx;

   std::cout << "Tol. Num Diffusion : " << tol_numdiffusion << std::endl;
   printf("Step Size (Delta t) :  %f\n", DELTAT);
   std::cout << "TSIZE : " << TSIZE << std::endl;
   std::cout << " Tsize-1*dt " << static_cast<double>(TSIZE-1)*DELTAT << " // Tsize*dt " << static_cast<double>(TSIZE)*DELTAT << std::endl;
   
   // CRIACAO DO MODELO
   Sedimentacao SedimentationPhysics(
       diameters, rho_p, rho_f, mu_f, settling_n, settling_lambda,
       Ne, tol_numdiffusion, phi_max     );

    /*
     * Class used to set BCs
     */
   std::vector<double> left_conditions;
   std::vector<double> right_conditions;
   left_conditions.assign(PDEDim, 0.0);
   right_conditions.assign(PDEDim, 0.0);
   
   std::vector<Bdry*> myBCs(2);
   myBCs[0] = new PrescribedFlux_L(PDEDim, CSTDim, SedimentationPhysics,
			left_conditions, dx);
		
            // PhysicsRCD& physicsRCD_p, std::vector<double> Feval_p, double dx_p,
            // double pDirichlet_p = 0.0, double uDirichlet_p = 0.0);

   myBCs[1] = new PrescribedFlux_R(PDEDim, CSTDim, SedimentationPhysics,
			right_conditions, dx);

   /*
     * Class for General discretization (In this case, many default values are
     *being used)
     */
   DiscretizationParallelized dct(
            SedimentationPhysics, // instância da classe physics
            myBCs,            // vetor de condições de fronteira
            PDEDim,           // quantidade de EDPs no sistema
            CSTDim,           // quantidade de restrições no sistema
            XSIZE,            // numero de pontos espaciais
            dx,               // distancia de um ponto consecutivo para outro
            false, // booleano que indica se há pressão (equação separada
                   // do sistema)
            false); // booleano que indica se há velocidade (equaçãoseparada do
                    // sistema)
                    
   PDEMemberConfig pdemc;
    pdemc.hasFlux() = true; // booleano que indica se a discretização
                            // considerará termo de fluxo
    pdemc.fluxIsUpwind()
            = false; // booleano que indica se o termo de fluxo é Upwind ou não
    pdemc.impParamFlux() = 0.5; // valor entre 0 e 1. É a constante que indica a
                                // média temporal usada para o termo de fluxo

    pdemc.hasDiffusion() = true;  // booleano que indica se a discretização
                                   // considerará a matrix de difusão

    pdemc.diffIsConstant() = true; // booleano que indica se a matriz de difusão
    // é composta por entradas constantes ou não
    pdemc.impParamDiffusion() = 0.5; // valor entre 0 e 1. É a constante que
                                     // indica a média temporal usada na
                                     // discretização do termo que contém a
                                     // matrix de difusão
    pdemc.hasReaction() = false;     // booleano que indica se a discretização
                                     // considerará termo de reação
    pdemc.impParamReaction() = 0.5;  // valor entre 0 e 1. É a constante que
                                     // indica a média temporal usada para o
                                     // termo de reação
    pdemc.impParamConstraint()
            = 1.0; // valor entre 0 e 1. É a constante que indica a média
    dct.setUpPDEMembers(pdemc); // temporal usada para o termo de restrição

    /*
     * Class for Newton Method
     */
    Newton nwt(dct,    // ponteiro para uma intância da classe Discretization
            NMAXITER); // valor máximo permitido de iterações para o método de
    // Newton

   /*
     * Preenchimento do vetor dct.x
     *
     * OBS: Isso até o momento deve ser feito pelo usuário
     */
    for (std::size_t i = 0; i != XSIZE; i++)
    {
        dct.x[i] = XLEFT + (static_cast<double>(i) ) *dx;
    }
    Set_idxcotas(dct.x, idx_cotas, cota_position, XSIZE, num_cotapos);
    for (int i=0; i< num_cotapos; i++)
		std::cout << "Cota ["<< i<<"] = "<< cota_position[i] << "   INDEX: "<< idx_cotas[i] 
		<< "  X: " << dct.x[idx_cotas[i]]<< std::endl;
    /*
     * The initial solution is calculated by this function
     */
    CalculateInitialSolution(dct.V, XSIZE, Ne, initial_concentrations);

	// CHECK CONSERVATION
	std::cout << "XSIZE in code: " << XSIZE << std::endl;
	double *initial_total_concentrations = new double[Ne];
	for (std::size_t i=0; i<Ne; i++)
	{
		std::cout << "---" << ( dx*static_cast<double>(XSIZE-1) ) << std::endl;
		initial_total_concentrations[i] = static_cast<double>(XSIZE-1) * initial_concentrations[i]*dx;
		std::cout << "Total Initial Particles Concentration " << i << " :: " 
				  << initial_total_concentrations[i] << std::endl;
	}
	double total_concentration;
    double total_ini_concentration; // Concentracao Inicial Total. Usado para calcular o topo.

	double DeltaT_Print;
	double nextPrint;
	double currentTime {0.0};
	bool found = false;
	size_t nprint {0};
	double sim_Deltat;
	
	DeltaT_Print =  TFINAL/static_cast<double>(NOUTPUTS-1);
	nextPrint = DeltaT_Print;
	if (DeltaT_Print < DELTAT )
		DeltaT_Print = DELTAT;
	std::cout << "DeltaT_Print = " << DeltaT_Print << "  DELTAT =" << DELTAT << std::endl;
	
	// SAVE INITIAL SOLUTION
	std::cout << " # Saving at Time: "<< currentTime  << "  nprint= "<< nprint << std::endl;
	Save_Solution(foldername, dct.x, dct.V, XSIZE, Ne, nprint, currentTime);
	nprint++;
	
    // CHECK IF TIME=0 IS A COTA_TIME AND SAVE
    std::string cotaname;
    cotaname = foldername + "/cotas.data";
    std::ofstream cotafile(cotaname); 
    cotafile.precision(16);
    cotafile << std::scientific;
    
    for (std::size_t i =0; i<num_cotapos ; i++)
        cotafile << dct.x[idx_cotas[i]] << " ";
    cotafile << std::endl;
    cotafile << "Cota Time |  Solution at Given Positions" << std::endl;
    cotafile << "0.000000 " ;
    for (std::size_t i=0; i<num_cotapos; i++)
    {
        double sumphi {0.0};
        for (int k = 0; k<Ne; k++)
            sumphi += dct.V[idx_cotas[i]*Ne + k] ;
        cotafile << sumphi << " ";
    }
    cotafile << std::endl;


    // SAVE SIMULATION DATA: 
    std::cout << "SAVING SIMULATION INFO" << std::endl;
    std::string simname;
    simname = foldername + "/simulation.info";
    std::cout << simname << std::endl;
    std::ofstream simufile(simname);
    simufile.precision(16);
    simufile << std::scientific;
    simufile << "Simulation Information" << std::endl;
    simufile << "Nx : " << XSIZE << std::endl;
    simufile << "Domain Length : " << LENGTH << std::endl;
    simufile << "Simulation Final Time : " << TFINAL << std::endl;
    simufile << "Delta T : " << DELTAT << std::endl;
    simufile << "**********************" << std::endl;
    simufile << "Number of Particles : " << Ne << std::endl;
    simufile << "Particle Num. | Initial Concentration | Diameter | Density "
             << std::endl;
    for (int i = 0; i < Ne; i++)
        simufile << initial_concentrations[i] << "  " << diameters[i] << "  "
                 << rho_p[i] << std::endl;
    simufile << "Maximal Total Concentration : " << phi_max << std::endl;
    simufile << "Hindered Settling Function - Parameter n : " << settling_n
             << std::endl;
    simufile << "Hindered Settling Function - Parameter lambda : "
             << settling_lambda << std::endl;
    simufile << "Fluid Viscosity : " << mu_f << std::endl;
    simufile << "Fluid Density   : " << rho_f << std::endl;
    simufile << "**********************************" << std::endl;
    simufile << "Number of Output Files in Time : " << NOUTPUTS << std::endl;
    simufile << "Number of Cota Positions to Save : " << num_cotapos
             << std::endl;
    simufile << "Cotas : " << std::endl;
    for (int i = 0; i < num_cotapos; i++)
        simufile << "      Cota x : " << cota_position[i] << std::endl;
    simufile << "**********************************" << std::endl;
    simufile << "Numerical Diffusion Factor : "
             << "TO BE PROPERLY DEFINED " << std::endl;
    simufile << "Numerical Diffusion Value: " << tol_numdiffusion << std::endl;
    simufile << "Newton - Max. Iterations : " << NMAXITER << std::endl;
    simufile << "Newton - Zero Tolerance  : " << ZEROTOL << std::endl;
    simufile << "Newton - DIFF TOL        : " << DIFFTOL << std::endl;
    simufile << "List of Saved Times" << std::endl;
    simufile << "0.000" << std::endl;

   
   /*
    * 
    * SIMULATION LOOP
    * 
    * */ 
   
   // Create file to save  the position of top of sedimentation ( t , x) 
   std::cout << "Abrindo arquivo para salvar o topo da sedimentacao" << std::endl;
   std::string topo_nome;
    topo_nome = foldername + "/topo.data";
    std::cout << topo_nome << std::endl;
    std::fstream arquivo_topo;
    arquivo_topo.open(topo_nome, std::fstream::out);

    arquivo_topo.precision(16);
    arquivo_topo << std::scientific; 

    //  tempo x posicao
    //arquivo_topo << currentTime << "  " << LENGTH << std::endl;
    scansave_topo(arquivo_topo, -1.0, dct.x, dct.V, XSIZE, Ne);



   int cota_print {0};
   for (std::size_t tc = 1; tc <= TSIZE; tc++) {
	   
	   if ( (currentTime+DELTAT)>TFINAL)
	   {
		   sim_Deltat =  TFINAL-currentTime;
	   }
	   else
	   {
		   sim_Deltat = DELTAT;
        //    if (current_cota<num_cotatimes) {
               
		//     if ( ((currentTime+DELTAT)>cota_time[current_cota]) && (currentTime<cota_time[current_cota]) )
        //     {
        //        if (fabs( (currentTime+DELTAT) - cota_time[current_cota])>1.0e-6)
        //        {
        //             std::cout<< "==================================" << std::endl;
        //             std::cout << "!!Adjust to Print Cota" << std::endl;
        //             cota_print = 1;
        //             sim_Deltat = cota_time[current_cota] - currentTime;
        //             std::cout << " New SimDeltat : " << sim_Deltat << std::endl;
        //             std::cout << "idx current cota: " << current_cota << std::endl;
        //             std::cout << "currentTime : " << currentTime << std::endl;
        //             std::cout << "currentTime+DeltaT : " << currentTime+DELTAT << std::endl;
        //             std::cout << "Cota Time: " << cota_time[current_cota] << std::endl;
        //             std::cout<< "==================================" << std::endl;
        //        }
        //     }
        //    }         
	   }

        if ((tc % 100) == 0)
        {
            //  std::cout << "\r> Calculating solution vector at time : "
            //            << currentTime+dt << std::flush;
            std::cout << "Computing sol. at time: " << currentTime+sim_Deltat << std::endl;            
	    }

        /*
         * The newton method is executed
         */
        switch (nwt.TimeStep(sim_Deltat, currentTime, ZEROTOL, DIFFTOL)) {
        case 0: {
            found = true;
            std::cerr << std::endl << "The process was killed because the linear system "
                         "was not solved."
                      << std::endl;
        } break;

        case -1: {
            found = true;
            std::cerr << std::endl << "The process was killed because the maximum number "
                         "of iterations was reached."
                      << std::endl;
        } break;
       
        } // switch end
        currentTime += sim_Deltat; // Update the time  
        

        //scansave_topo(arquivo_topo, currentTime, dct.x, dct.V, XSIZE, Ne);


        // PLOT AT GIVEN COTA TIME...
        // if (current_cota<num_cotatimes)
        // {
        //     if ((cota_print) || (fabs(currentTime - cota_time[current_cota])<1.0e-6))
        //     {
        //     std::cout << "Saving At Cota " << currentTime << "  -  "<< cota_time[current_cota] << std::endl;
        //     //std::cout << "   1) " << (fabs(currentTime - cota_time[current_cota]) )<< std::endl;
        //     //std::cout << cota_print << std::endl;
        //     cotafile << cota_time[current_cota] << " ";
        //     for (std::size_t i =0; i<num_cotapos ; i++) 
        //     {
        //         double sumphi {0.0};
        //         for (int k = 0; k<Ne; k++)
        //             sumphi += dct.V[idx_cotas[i]*Ne + k] ;
        //         cotafile << sumphi << " ";
        //     }    
        //     cotafile << std::endl;
            
        //     cota_print=0;
        //     current_cota++;
        //     }
        // }

		// SAVE TIME STEP SOLUTION + cota:
		if ( (currentTime > nextPrint) && (nprint< (NOUTPUTS-1)) )
		{
			std::cout << " # Saving at Time: "<< currentTime  << "  nprint= "<< nprint << std::endl;
			std::cout << " @@@ Save Solution File ..." << std::endl;
            Save_Solution(foldername, dct.x, dct.V, XSIZE, Ne, nprint, currentTime);
            
            std::cout << " @@@ Save Solution Time to file simufile..." << std::endl;
            simufile << currentTime << std::endl;

            std::cout << " @@@ Salvando Topo da Sedimentacao no Arquivo" << std::endl;
            scansave_topo(arquivo_topo, currentTime, dct.x, dct.V, XSIZE, Ne);
             

            std::cout << " ##### Saving at Cota" << std::endl;
            cotafile << currentTime << " ";
            for (std::size_t i =0; i<num_cotapos ; i++) 
            {
                double sumphi {0.0};
                for (int k = 0; k<Ne; k++)
                    sumphi += dct.V[idx_cotas[i]*Ne + k] ;
                cotafile << sumphi << " ";
            } 
            cotafile << std::endl;
			nextPrint += DeltaT_Print;
			nprint++;
		}
		
		// SAVE SOLUTION OF X POSITION


		// CONSERVARTION OF MASS CHECKING
		// if ((tc % 1000) == 0)
		// {
		// 	std::cout << "Checking Conservation: " << std::endl;
		// 	for (std::size_t i=0; i<Ne; i++)
		// 	{
		// 		total_concentration = 0.0;
		// 		for (std::size_t k=0; k<(XSIZE-1); k++)
		// 		{
		// 			total_concentration += ( dct.V[k*Ne + i] + dct.V[(k+1)*Ne+i] )*0.5*dx;
		// 		}
		// 		std::cout << "Total of Particle [ "<<i<<" ]  : " << total_concentration  
		// 		<< " ||  Erro Relativo e Absoluto:  " << 
		// 		(total_concentration - initial_total_concentrations[i])/total_concentration 
		// 		<< " %  &  " << (total_concentration - initial_total_concentrations[i]) << std::endl;
		// 	}
		// }

        if (found)
            break;
    }

	std::cout << " # Saving at Time: "<< currentTime  << "  nprint= "<< nprint << std::endl;
	Save_Solution(foldername, dct.x, dct.V, XSIZE, Ne, nprint, currentTime);
    simufile << currentTime << std::endl;
    std::cout << " ##### Saving at Cota" << std::endl;
    
    cotafile << currentTime << " ";
    for (std::size_t i =0; i<num_cotapos ; i++) 
    {
        double sumphi {0.0};
        for (int k = 0; k<Ne; k++)
            sumphi += dct.V[idx_cotas[i]*Ne + k] ;
        cotafile << sumphi << " ";
    } 
    cotafile << std::endl;
    cotafile.close();

    simufile.close();
    arquivo_topo.close();

    delete [] cota_time;

   delete myBCs[0];
   delete myBCs[1];   
   
   return 0;
}

