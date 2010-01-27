inline double fround(double n, unsigned d)
{
  return floor(n * pow(10., d) + .5) / pow(10., d);
}
