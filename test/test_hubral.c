#include "Unity/unity.h"
#include "forward.h"
#include "rungekutta.h" 
#include <stdio.h>
#include <rsf.h>

int n[2]={301,1001};
float d[2]={0.01,0.01};
float o[2]={0.,-2.};
float v0=1.508;
float *slow;
float *slow2;

void init(){
	int i, j;
	float z;
	float ss=1./(v0*v0);
	slow = sf_floatalloc(n[0]*n[1]);

	for(j=0; j<n[1]; j++){
		for(i=0; i<n[0]; i++){
			z=i*d[0];
			if(z<1.){
				slow[i+j*n[0]] = 1.508;
			}else{
				slow[i+j*n[0]] = 2.0;
			}
		}
	}
	for(i=0;i<n[0]*n[1];i++){
		slow[i]=1./(slow[i]*slow[i]);
	}

}

void init2()
/*< Stack layers model >*/
{
	int i, j;
	float z;
	float ss=1./(v0*v0);
	slow2 = sf_floatalloc(n[0]*n[1]);

	for(j=0; j<n[1]; j++){
		for(i=0; i<n[0]; i++){
			z=i*d[0];
			if(z<1.){
				slow2[i+j*n[0]] = 1.508;
			}else if(z<1.85){
				slow2[i+j*n[0]] = 1.69;
			}else if(z<2.5){
				slow2[i+j*n[0]] = 1.75;
			}else{
				slow2[i+j*n[0]] = 2.0;
			}
		}
	}
	for(i=0;i<n[0]*n[1];i++){
		slow2[i]=1./(slow2[i]*slow2[i]);
	}

}


void setUp(){}

void tearDown(){}

void test_getRNIPUsingHubralLawsInConstantVelocityModel()
/*< For a constant velocity layer without interfaces, RNIP is simply equal to normal ray length
(propagation law). So, this test propagates normal rays in the first layer only and calculates RNIP
using Hubral's propagation law.
Example: For a 1km length layer, RNIP should be equal to 1.
>*/
{
        float x[2];
        float p[2];
        raytrace rt;
        float **traj;
        float normalRayAngleRad;
	float a[1]={0.,};
        int nt=10000;
        float dt=0.001;
        int it;
        float t;
        int is;
	float **s;
	float rnip=0.;
	int i;
	float sz[5]={1.,1.,1.,1.,1.};
	int nsz=5;
	float osz=0., dsz=1.;
	float vv[2]={1.508,2.};

	s = sf_floatalloc2(2,1);

        /* initialize ray tracing object */
        rt = raytrace_init(2,true,nt,dt,n,o,d,slow,ORDER);
        traj = sf_floatalloc2(2,nt+1);
	normalRayAngleRad = a[0]*DEG2RAD;

	for(i=0;i<5;i++){
		s[0][0] = 1.;
		s[0][1] = i*0.5+2.5;
		/* Set initial ray point and ray vector */
		setInitialRayPointAndRayVector(s,x,p,0,normalRayAngleRad);

		/* Ray tracing */
		it = trace_ray (rt, x, p, traj);

		rnip = calculateRNIPWithHubralLaws(rt,traj,it,vv,2,0,sz,nsz,osz,dsz);

		TEST_ASSERT_FLOAT_WITHIN(0.1,1.0,rnip);
	}

	raytrace_close(rt);
        free(traj);
	free(s);
}

void test_getRNIPUsingHubralLawsInConstantVelocityModelSecondLayer()
/*< REPEAT THE FIRST TEST FOR THE SECOND LAYER
For a constant velocity layer without interfaces, RNIP is simply equal to normal ray length
(propagation law). So, this test propagates normal rays in the first layer only and calculates RNIP
using Hubral's propagation law.
Example: For a 2km length layer, RNIP should be equal to 2.
>*/
{
        float x[2];
        float p[2];
        raytrace rt;
        float **traj;
        float normalRayAngleRad;
	float a[1]={0.,};
        int nt=1000;
        float dt=0.001;
        int it;
        float t;
        int is;
	float **s;
	float rnip=0.;
	int i;
	float vv[2]={1.508,2.};
	float sz[1]={1.};
	int nsz=1;
	float osz=0., dsz=1.;

	s = sf_floatalloc2(2,1);

        /* initialize ray tracing object */
        rt = raytrace_init(2,true,nt,dt,n,o,d,slow,ORDER);
        traj = sf_floatalloc2(2,nt+1);
	normalRayAngleRad = a[0]*DEG2RAD;

	for(i=0;i<5;i++){
		s[0][0] = 3.;
		s[0][1] = i*0.5+2.5;
		/* Set initial ray point and ray vector */
		setInitialRayPointAndRayVector(s,x,p,0,normalRayAngleRad);

		/* Ray tracing */
		it = trace_ray (rt, x, p, traj);

		rnip = calculateRNIPWithHubralLaws(rt,traj,nt,vv,2,1,sz,nsz,osz,dsz);

		TEST_ASSERT_FLOAT_WITHIN(0.1,2.,rnip);
	}

	raytrace_close(rt);
        free(traj);
	free(s);
}

void test_getRNIPUsingHubralLawsInTwoLayersVelocityModel()
/*< REPEAT THE FIRST TEST FOR THE SECOND LAYER
For a constant velocity layer without interfaces, RNIP is simply equal to normal ray length
(propagation law). So, this test propagates normal rays in the first layer only and calculates RNIP
using Hubral's propagation law.
Example: For a 2km length layer, RNIP should be equal to 2.
>*/
{
        float x[2];
        float p[2];
        raytrace rt;
        float **traj;
        float normalRayAngleRad;
	float a[1]={0.,};
        int nt=5000;
        float dt=0.001;
        int it;
        float t;
        int is;
	float **s;
	float rnip=0.;
	int i;
	float vv[2]={1.508,2.};
	float sz[5]={1.,1.85,1.,1.,1.};
	int nsz=5;
	float osz=0., dsz=1.;

	s = sf_floatalloc2(2,1);

        /* initialize ray tracing object */
        rt = raytrace_init(2,true,nt,dt,n,o,d,slow,ORDER);
        traj = sf_floatalloc2(2,nt+1);
	normalRayAngleRad = a[0]*DEG2RAD;

	for(i=0;i<1;i++){
		s[0][0] = 3.;
		s[0][1] = i*0.5+2.5;
		/* Set initial ray point and ray vector */
		setInitialRayPointAndRayVector(s,x,p,0,normalRayAngleRad);

		/* Ray tracing */
		it = trace_ray (rt, x, p, traj);

		rnip = calculateRNIPWithHubralLaws(rt,traj,it,vv,2,1,sz,nsz,osz,dsz);
		//sf_warning("%d rnip=%f %f %f",it,rnip,traj[it-1][1],traj[it-1][0]);

		TEST_ASSERT_FLOAT_WITHIN(0.1,3.65,rnip);
	}

	raytrace_close(rt);
        free(traj);
	free(s);
}

void test_getRNIPUsingHubralLawsInStackLayersVelocityModel()
/*< REPEAT THE FIRST TEST FOR THE SECOND LAYER
For a constant velocity layer without interfaces, RNIP is simply equal to normal ray length
(propagation law). So, this test propagates normal rays in the first layer only and calculates RNIP
using Hubral's propagation law.
Example: For a 2km length layer, RNIP should be equal to 2.
>*/
{
        float x[2];
        float p[2];
        raytrace rt;
        float **traj;
        float normalRayAngleRad;
	float a[1]={0.,};
        int nt=5000;
        float dt=0.001;
        int it;
        float t;
        int is;
	float **s;
	float rnip=0.;
	int i;
	float vv[4]={1.508,1.69,1.75,2.};
	float sz[4]={1.,1.85,2.5,3.};
	int nsz=5;
	float osz=0., dsz=1.;

	s = sf_floatalloc2(2,1);

        /* initialize ray tracing object */
        rt = raytrace_init(2,true,nt,dt,n,o,d,slow2,ORDER);
        traj = sf_floatalloc2(2,nt+1);
	normalRayAngleRad = a[0]*DEG2RAD;

	// First interface
	for(i=0;i<1;i++){
		s[0][0] = 1.;
		s[0][1] = i*0.5+2.5;
		/* Set initial ray point and ray vector */
		setInitialRayPointAndRayVector(s,x,p,0,normalRayAngleRad);

		/* Ray tracing */
		it = trace_ray (rt, x, p, traj);

		rnip = calculateRNIPWithHubralLaws(rt,traj,it,vv,4,0,sz,nsz,osz,dsz);
		//sf_warning("%d rnip=%f %f %f",it,rnip,traj[it-1][1],traj[it-1][0]);

		TEST_ASSERT_FLOAT_WITHIN(0.1,1.,rnip);
	}

	// Second interface
	for(i=0;i<1;i++){
		s[0][0] = 1.85;
		s[0][1] = i*0.5+2.5;
		/* Set initial ray point and ray vector */
		setInitialRayPointAndRayVector(s,x,p,0,normalRayAngleRad);

		/* Ray tracing */
		it = trace_ray (rt, x, p, traj);

		rnip = calculateRNIPWithHubralLaws(rt,traj,it,vv,4,1,sz,nsz,osz,dsz);
		//sf_warning("%d rnip=%f %f %f",it,rnip,traj[it-1][1],traj[it-1][0]);

		TEST_ASSERT_FLOAT_WITHIN(0.1,1.95,rnip);
	}

	// Third interface
	for(i=0;i<1;i++){
		s[0][0] = 2.5;
		s[0][1] = i*0.5+2.5;
		/* Set initial ray point and ray vector */
		setInitialRayPointAndRayVector(s,x,p,0,normalRayAngleRad);

		/* Ray tracing */
		it = trace_ray (rt, x, p, traj);

		rnip = calculateRNIPWithHubralLaws(rt,traj,it,vv,4,2,sz,nsz,osz,dsz);
		//sf_warning("%d rnip=%f %f %f",it,rnip,traj[it-1][1],traj[it-1][0]);

		TEST_ASSERT_FLOAT_WITHIN(0.1,2.7,rnip);
	}

	// Fourth interface
	for(i=0;i<1;i++){
		s[0][0] = 3.;
		s[0][1] = i*0.5+2.5;
		/* Set initial ray point and ray vector */
		setInitialRayPointAndRayVector(s,x,p,0,normalRayAngleRad);

		/* Ray tracing */
		it = trace_ray (rt, x, p, traj);

		rnip = calculateRNIPWithHubralLaws(rt,traj,it,vv,4,3,sz,nsz,osz,dsz);

		TEST_ASSERT_FLOAT_WITHIN(0.1,3.36,rnip);
	}


	raytrace_close(rt);
        free(traj);
	free(s);
}

int main(int argc, char* argv[]){

	init();
	init2();
        UNITY_BEGIN();
        RUN_TEST(test_getRNIPUsingHubralLawsInConstantVelocityModel);
        RUN_TEST(test_getRNIPUsingHubralLawsInConstantVelocityModelSecondLayer);
        RUN_TEST(test_getRNIPUsingHubralLawsInTwoLayersVelocityModel);
	RUN_TEST(test_getRNIPUsingHubralLawsInStackLayersVelocityModel);
        return UNITY_END();
}

