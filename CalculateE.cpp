#include <iostream>
using namespace std;

int main()
{
  
  double a, b, cH, E;

  cout << "Give a and b for E = a*# + b" << endl;
  cin >> a >> b;

  cout << "Give Energy in channels: " << endl;

  cin >> cH;
  
  while(!cH==0)
  {
    E = a*cH+b;
    printf("Channel : %.10f  -->  Energy: %.1f keV \n",cH,E);     
    cout << "Give next one, type '0' to quit: " << endl;
    cin >> cH;
  }
 

return 1;

}
