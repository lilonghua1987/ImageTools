#include <iostream>
#include <cstdlib>
#include <cmath>
#include "MyPSO.h"
using namespace std;

Random Particle::_rnd;
Random ParticleSwarm::_rnd;

int main()
{
	cout<<"PSO Begaining..."<<endl;

	int iteration = 0;
	int countSame = 0;
	double minimum = 1.0e+100;
	double oldMinimum = 1.0e+100;

	SphereFunction fuction;
	int swarmSize = 100;
	int dimension = 10;

	FunctionMinimizingParticleSwarm swarm(&fuction,swarmSize,dimension);

	while (countSame < 25)
	{
		swarm.Iteration();
		minimum = swarm.BestCost;
		iteration++;

		if (abs(minimum - oldMinimum) < 0.00001)
		{
			countSame++;
		}
		else
		{
			countSame = 0;
		}
		oldMinimum = minimum;

		cout<<iteration<<"\t"<<swarm.CurrentBestCost()<<"\n"<<endl;

	}

	cout<<"Best Cost Are:"<<swarm.BestCost<<"\n"<<endl;
	cout<<"Best Position Are:\n"<<endl;
	vector<double>::iterator iter;
	for (iter = swarm.BestPosition.begin(); iter != swarm.BestPosition.end(); ++iter)
	{
		cout<<*iter<<" "<<endl;
	}

	system("pause");
	return 1;
}