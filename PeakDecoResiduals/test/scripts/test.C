void setArray(vector<float> &x){
  x.clear();
  x.push_back(1.3);
  x.push_back(2.4);
}

int run(){

  vector<float> y;
  setArray(y);
  cout<<y.at(0)<<" "<<y.at(1)<<endl;


//   cout<<"Start"<<endl;
//   float y[]=getArray();
//   cout<<"Got array"<<endl;
//   cout<<y[0]<<" "<<y[1]<<" "<<y[2]<<endl;
//   cout<<"Finished"<<endl;
//   delete y[];

  return 0;
}

// float *getArray(){
//   float x[]={1.,14,34.54};
//   return x;
// }

// int run(){
//   cout<<"Start"<<endl;
//   float *y=getArray();
//   cout<<"Got array"<<endl;
//   cout<<y[0]<<" "<<y[1]<<" "<<y[2]<<endl;
//   cout<<"Finished"<<endl;


//   return 0;
// }



/*
#include <iostream>

int *foo(int *array)

      {

      int *start = array;

      while ( *array )

      {

      std::cout << *array++ << ' ';

      }

      std::cout << std::endl;

      return start;

      }

       

      int *bar(int *array)

      {
      int *start = array;
  17.
      while ( *array )
  18.
      {
  19.
      *array++ += 1;
  20.
      }
  21.
      return start;
  22.
      }
  23.
       
  24.
      int main()
  25.
      {
  26.
      int myarray[] = {1,2,3,4,5,0};
  27.
      foo(bar(foo(bar(foo(myarray)))));
  28.
      return 0;
  29.
      }

*/
