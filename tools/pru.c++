#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

#define KORD 5
#define L_ITERVALS 20

using namespace std;

int main(void){

	unsigned int nb = KORD+L_ITERVALS-3;

	ifstream in1("salida1");
	ifstream in2("salida2");
	double d1, d2;
	int i=0;
	while(in1>>d1 &&  in2>>d2)
	{
		//double d1 = std::stod(s1);
		//double d2 = std::stod(s2);
		
		if(fabs(d1 - d2) > 1e-5){
			cout << i << " " << d1 << " " << d2 << endl ;
		}
		i++;
	}
	if(in2>>d2){
		cout << "el segundo era mayor que el primero" << endl;
	}
	else if(in1 >> d1){
		cout << " el primero era mayor que el segundo" << endl;
	}
	in1.close();
	in2.close();
	return 0;
}