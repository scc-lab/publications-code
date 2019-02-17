
// ========================================================================
//  QMotor - A PC Based Real-Time Graphical Control Environment
//  (c) 2000 QRTS
//  ---------------------------------------------------------------
//  Control Program : SATRISE.cpp
//  Description     : Saturated RISE controller
//  Author          : Nic Fischer
//  Last Modified   : 2011.
//  For more information on this code, please see the paper:
//  N. Fischer, Z. Kan, W. E. Dixon, "Saturated RISE Feedback Control for Euler-Lagrange Systems,"
//	American Control Conference, Ontario, CA, 2012.
// ========================================================================

// ----- QRTS libraries -----
#include "ControlProgram.hpp"
#include "IOBoardClient.hpp"
#include "ButterworthFilter.hpp"

// ----- C standard libraries -----
#include <math.h>
#include <iostream.h>
#include <stdio.h>
#include <iomanip.h>
#include <stdlib.h>

# define pi 3.141592654
# define n 2 // number of links in the robot
//============================================================================
// Class definition of the SATRISE class
//============================================================================

class SATRISE : public ControlProgram
{
	protected:

        // Error system variables
		double *q;
		double *q_dot;
		double *qd;
		double *qd_dot;
		double *e1;
		double *e1_dot;
		double *e2;
		double *ef;
		double *ef_dot;
        double *v;
        double *v_dot;
		double *torque;
		double *voltage;
		double *q_pulse;

        // Control gains
		double *gamma1;
		double *gamma2;
		double *alpha1;
		double *alpha2;
		double *alpha3;
		double *beta;
		double threshold;

        // Integration helper terms
		double *e1_prev;
		double *q_prev;
		double *q_dot_prev;
        double *ef_prev;
		double *ef_dot_prev;
		double *v_prev;
		double *v_dot_prev;

		// Saturation terms
		double sat_ef;
		double sat_v;

		// Other random variables
		double t;
        double M_denom;
        double *Minv;
        double *p;
        FILE * pFile;

        // Output variables
        double *e1_deg;
        double *q_deg;
        double *qd_deg;


		ButterworthFilter<double> *filter;

		// Clients
		IOBoardClient *iobc1, *iobc2;

    public:

        SATRISE(int argc, char *argv[]) : ControlProgram (argc, argv){};
        ~SATRISE () {};

        void allocateVariables();
        void registerLogVariables();
        void registerControlParameters();
        void deleteVariables();
        void initializeVariables();
        int signum(double x);

        // Default QMotor override functions
        virtual int enterControl();
		virtual int startControl();
		virtual int control();
		virtual int stopControl();
		virtual int exitControl();

};

void SATRISE::allocateVariables()
{
    // Allocate all variables in the class
    q = new double[2];
    q_dot = new double[2];
    qd = new double[2];
    qd_dot = new double[2];
    e1 = new double[2];
    e1_dot = new double[2];
    e2 = new double[2];
    ef = new double[2];
    ef_dot = new double[2];
    v = new double[2];
    v_dot = new double[2];
    torque = new double[2];
    voltage = new double[2];
    q_pulse = new double[2];

    gamma1 = new double[2];
    gamma2 = new double[2];
    alpha1 = new double[2];
    alpha2 = new double[2];
    alpha3 = new double[2];
    beta = new double[2];

    e1_prev = new double[2];
    q_prev = new double[2];
    q_dot_prev = new double[2];
    ef_prev = new double[2];
    ef_dot_prev = new double[2];
    v_prev = new double[2];
    v_dot_prev = new double[2];

    Minv = new double[3];
    p = new double[3];

    e1_deg = new double[2];
    q_deg = new double[2];
    qd_deg = new double[2];

    filter = new ButterworthFilter<double> [2];
}

void SATRISE::deleteVariables()
{
    delete [] q;
    delete [] q_dot;
    delete [] qd;
    delete [] qd_dot;
    delete [] e1;
    delete [] e1_dot;
    delete [] e2;
    delete [] ef;
    delete [] ef_dot;
    delete [] v;
    delete [] v_dot;
    delete [] torque;
    delete [] voltage;
    delete [] q_pulse;

    delete [] gamma1;
    delete [] gamma2;
    delete [] alpha1;
    delete [] alpha2;
    delete [] alpha3;
    delete [] beta;

    delete [] e1_prev;
    delete [] q_prev;
    delete [] q_dot_prev;
    delete [] ef_prev;
    delete [] ef_dot_prev;
    delete [] v_prev;
    delete [] v_dot_prev;

    delete [] Minv;
    delete [] p;

    delete [] e1_deg;
    delete [] q_deg;
    delete [] qd_deg;

    delete [] filter;
}

void SATRISE::registerLogVariables()
{
    // Register all variables that we want to log
    registerLogVariable(q_deg, "q", "Current Position (Deg)", 2);
    registerLogVariable(qd_deg, "qd", "Desired Position (Deg)", 2);
	registerLogVariable(e1_deg, "e1", "Error 1", 2);
	registerLogVariable(e2, "e2 (rad/s)", "Error 2", 2);
	registerLogVariable(torque, "torque", "Torque", 2);
	registerLogVariable(v, "v", "Rise Integral", 2);
	registerLogVariable(ef, "ef", "Filtered error", 2);
	registerLogVariable(q_dot, "q_dot", "Current Velocity", 2);
	registerLogVariable(voltage, "voltage", "Voltage", 2);
}

void SATRISE::registerControlParameters()
{
    registerControlParameter (gamma1, "gamma1", "gamma1", 2);
    registerControlParameter (gamma2, "gamma2", "gamma2", 2);
    registerControlParameter (alpha1, "alpha1", "alpha1", 2);
    registerControlParameter (alpha2, "alpha2", "alpha2", 2);
    registerControlParameter (alpha3, "alpha3", "alpha3", 2);
    registerControlParameter (beta, "beta", "beta", 2);
    registerControlParameter (&threshold, "threshold", "Tolerance");
    registerControlParameter (&sat_ef, "sat_ef", "Saturation value for ef", 2);
    registerControlParameter (&sat_v, "sat_v", "Saturation value for v", 2);

}

void SATRISE::initializeVariables()
{
    ef[0] = 0; ef[1] = 0;
    v[0] = 0; v[1] = 0;

    e1_prev[0] = 0; e1_prev[1] = 0;
    q_prev[0] = 0; q_prev[1] = 0;
    q_dot_prev[0] = 0; q_dot_prev[1] = 0;
    ef_prev[0] = 0; ef_prev[1] = 0;
    ef_dot_prev[0] = 0; ef_dot_prev[1] = 0;
    v_prev[0] = 0; v_prev[1] = 0;
    v_dot_prev[0] = 0; v_dot_prev[1] = 0;

    p[0]=5.473; p[1]=2.193; p[2]=1.242;

    t = 0;

    pFile=fopen("log.txt","w");

}

int SATRISE::signum(double x)
{
	return x<0 ? -1 : x>0;
}


//============================================================================
// SATRISE::enterControl
// --------------------------------------------------------------------------
// This function is called when the control program is loaded. In standalone
// mode, this happens immediately. When using the GUI, it happens when the
// user loads the control program.
//============================================================================
int SATRISE::enterControl()
{
    allocateVariables();
    registerLogVariables();
	registerControlParameters();

	// Set all control parameters initially to zero
	//clearAllControlParameters();

	// Start message
	d_messageStream
		<< endl << "----- " << d_applicationName << " -----" << endl
		<< "This is SATRISE controller for a Two-Link robot. Developed by Nic Fischer, 2011." << endl << endl;

	return 0;
}


//============================================================================
// SATRISE::startControl
// --------------------------------------------------------------------------
// Called each time a control run is started. If running from the GUI, this
// will be called each time the START button is pushed.
//============================================================================

int SATRISE::startControl()
{
	initializeVariables();
	clearAllLogVariables();

	// ----- Initialize your clients here	-----
	char *ioboardServerName1 = "qrts/iobs0";
	iobc1 = new IOBoardClient(ioboardServerName1);

	// Reset encoders
	iobc1->setEncoderValue(0,0);

	if (iobc1->isStatusError())
	{
		d_status.setStatusError()
			<< d_applicationName << ": [SATRISE::startControl()] "
			<< "Error connecting to IO board server " << ioboardServerName1 << endl;
		delete iobc1;
		iobc1 = 0;
		return -1;
	}

	char *ioboardServerName2 = "qrts/iobs0";   //CHECK IF IT HAS TO BE IOBS0 OR IOBS1?
	iobc2 = new IOBoardClient(ioboardServerName2);

	// Reset encoders
	iobc2->setEncoderValue(1,0);

	if (iobc2->isStatusError())
	{
		d_status.setStatusError()
			<< d_applicationName << ": [SATRISE::startControl()] "
			<< "Error connecting to IO board server " << ioboardServerName2 << endl;
		delete iobc2;
		iobc2 = 0;
		return -1;
	}

    for(int i=0; i<2; i++)
    {
        filter[i].setCutOffFrequency(100000);
        filter[i].setSamplingTime(d_controlPeriod);
        filter[i].setDampingRatio(1.0);
        filter[i].setAutoInit();
    }

    return 0;

}


//============================================================================
// SATRISE::control
//----------------------------------------------------------------------------
// Called each control cycle. Do your input, control computations, and output
// here. If you return 0, the control will continue to execute. If you return
// nonzero, the control will abort. You may want to abort if some error
// condition occurs (excessive velocity,  etc.)
//============================================================================

int SATRISE::control()
{
    // Time variable
    t = d_elapsedTime;

	// Assumptions:
	// Inital velocities are assumed to be zero.
	// Initial position may or may not be zero.
	// In this case, initial desired position is zero and initial actual position is not zero.

	//Desired Trajectory in radians
   	qd[0] = 60*(pi/180)*sin(2*t)*(1-exp(-0.01*pow(t,3)));
   	qd[1] = 60*(pi/180)*sin(2*t)*(1-exp(-0.01*pow(t,3)));
   	qd_dot[0] = (6283*pow(t,2)*sin(2*t))/(200000*exp(pow(t,3)/100)) - (6283*cos(2*t)*(1/exp(pow(t,3)/100) - 1))/3000;
   	qd_dot[1] = (6283*pow(t,2)*sin(2*t))/(200000*exp(pow(t,3)/100)) - (6283*cos(2*t)*(1/exp(pow(t,3)/100) - 1))/3000;
	if(t>15)
	{
		// Go home!
		for(int i=0; i<2; i++)
		{
			qd[i] = qd[i]*exp(-(t-15));
			qd_dot[i] = qd_dot[i]*exp(-(t-15));
		}
	}

	// Feedback from the motor (no of revolutions)
	q_pulse[0] = iobc1->getEncoderValue(0)/614400.0; //(614400 pulses/revolution)
	q_pulse[1] = iobc2->getEncoderValue(1)/614400.0; //(614400 pulses/revolution)

	// Actual trajectory in radians
	for(int i=0; i<2; i++)
	{
        q[i] = q_pulse[i]*2*pi;
	}

//	if (d_elapsedTicks !=0)
//	{
//	initq = 0;
//	}

    // Calculate M
    M_denom = pow(p[2],2)*cos(2*q[1]) - 2*p[0]*p[1] + 2*pow(p[1],2) + pow(p[2],2);
    Minv[0] = -(2*p[1])/M_denom;
    Minv[1] = (2*p[1]+2*p[2]*cos(q[1]))/M_denom;
    Minv[2] = (2*p[1]+2*p[2]*cos(q[1]))/M_denom;


    for(int i=0; i<2; i++)
    {
        //Find derivative
        q_dot[i] = filter[i].filter((q[i]-q_prev[i])/d_controlPeriod);
        q_prev[i] = q[i];

        // Calculate error terms
        e1[i] = qd[i]-q[i];
    }
        e2[0] = (qd_dot[0]-q_dot[0]) + alpha1[0]*tanh(e1[0]) + Minv[0]*tanh(e1[0]) + Minv[1]*tanh(e1[1]);
        e2[1] = (qd_dot[1]-q_dot[1]) + alpha1[1]*tanh(e1[1]) + Minv[1]*tanh(e1[0]) + Minv[2]*tanh(e1[1]);
        ef_dot[0] = pow(cosh(ef[0]),2)*(-gamma1[0]*e2[0] + Minv[0]*tanh(e1[0]) + Minv[1]*tanh(e1[1]) - gamma2[0]*tanh(ef[0]));
        ef_dot[1] = pow(cosh(ef[1]),2)*(-gamma1[1]*e2[1] + Minv[1]*tanh(e1[0]) + Minv[2]*tanh(e1[1]) - gamma2[1]*tanh(ef[1]));

    for(int i=0; i<2; i++)
    {
        //Aux control terms
        v_dot[i] = pow(cosh(v[i]),2)*(alpha2[i]*tanh(e2[i]) + beta[i]*signum(e2[i]) + alpha3[i]*e2[i] + alpha1[i]*pow(cosh(e1[i]),-2)*e2[i] - /*Mterm*/ + gamma2[i]*e2[i]);

    		// Normal RISE test
	   	//e2[i] = (qd_dot[i]-q_dot[i])+alpha1[i]*e1[i];
        	//v_dot[i] = (gamma1[i]+1)*alpha2[i]*e2[i] + beta[i]*signum(e2[i]);

        // Integrate aux terms
        ef[i] = ef_prev[i] + 0.5*d_controlPeriod*(ef_dot[i]+ef_dot_prev[i]);
        v[i] = v_prev[i] + 0.5*d_controlPeriod*(v_dot[i]+v_dot_prev[i]);

        // Check for saturation of aux terms
        if(ef[i]>sat_ef)
            ef[i] = sat_ef;
        if(ef[i]<-sat_ef)
            ef[i] = -sat_ef;
        if(v[i]>sat_v)
            v[i] = sat_v;
        if(v[i]<-sat_v)
            v[i] = -sat_v;
    }



    for (int i=0; i<2; i++)
    {
        // Calulate control
        	torque[i] = gamma1[i]*tanh(v[i]);
        	// Normal RISE
		//torque[i] = (gamma1[i]+1)*e2[i] + v[i];

        // Save previous values
        ef_prev[i] = ef[i];
        ef_dot_prev[i] = ef_dot[i];
        v_prev[i] = v[i];
        v_dot_prev[i] = v_dot[i];
    }

	//Convert torque to voltage
	voltage[0] = torque[0]*10/240;
	voltage[1] = torque[1]*10/20;

	// Sending the output voltage to DAC
	iobc1->setDacValue(0, voltage[0]);
	iobc2->setDacValue(1, voltage[1]);

	// Convert output variables from radians to degrees for plotting
	for(int i=0; i<2; i++)
	{
		e1_deg[i] = e1[i]*180/pi;
		q_deg[i] = q[i]*180/pi;
		qd_deg[i] = qd[i]*180/pi;
	}

	fprintf(pFile,"t: %f\t e1: %f %f\t e2: %f %f v_dot: %f %f\t ef_dot: %f %f\n", t,e1[0],e1[1],e2[0],e2[1],v_dot[0],v_dot[1],ef_dot[0],ef_dot[1]);

	return 0;
}

//============================================================================
// SATRISE::stopControl()
//----------------------------------------------------------------------------
// Called each time a control run ends. If running from the GUI, this
// will be called each time the STOP button is pushed, or when a timed run
// ends, or when the control aborts itself.
//============================================================================

int SATRISE::stopControl()
{
	// Zero out the DAC
	iobc1->setDacValue(0, 0);
	iobc2->setDacValue(1, 0);

	delete iobc1;
	delete iobc2;

	fclose(pFile);

	return 0;
}


//============================================================================
// SATRISE::exitControl
// --------------------------------------------------------------------------
// This function is called when the control is unloaded. In standalone
// mode, this happens after one control run has completed. When using
// the GUI, it happens when the user loads a new control program or
// exits the GUI.
//============================================================================

int SATRISE::exitControl()
{
	return 0;
}

//============================================================================
// main()
//----------------------------------------------------------------------------
// The main function instantiates the object and goes into the mainloop
//============================================================================

using namespace std; //from Vilas's code, i don't know what it is?

main (int argc, char *argv[])
{
	SATRISE cp(argc, argv);
	cp.run();
}

