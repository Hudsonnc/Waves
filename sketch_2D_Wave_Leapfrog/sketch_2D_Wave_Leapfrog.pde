int N = 70;
double t = 0;
double dx = 1/(double)N;
double dy = 1/(double)N;
double dt;
double c = 1;
double b = 0;//.001;
double[][] vals = new double[N][N];
double[][] newvals = new double[N][N];
double[][] oldvals = new double[N][N];

double pressure(double t) {
  return 0; //10000*cos(2*PI*(float)t*440);
}

void deepcopy(double[][] src, double[][] dest) {
  for (int i = 0; i < src.length; i++) {
    arrayCopy(src[i], dest[i]);
  }
}
 
void setup() {
  dt = .1*1/c * 1/sqrt(1/sq((float)dx) + 1/sq((float)dy));
  println(dt);
  double Cx = c*dt/dx;
  double Cy = c*dt/dy;
  println(sq((float)Cx)+sq((float)Cy));
  frameRate(240);
  print(frameRate);
  size(1000,1000, P3D);
  smooth();
  for (int i=0; i<N; i++) {
    for (int j=0; j<N; j++) {
      vals[i][j] = .5*exp(-.02*(pow((i-3*N/4),2) + pow((j-N/2),2))) + .5*exp(-.02*(pow((i-N/4),2) + pow((j-N/2),2)));
    }
  }
  deepcopy(vals, oldvals); //Enforce derivative at t=0 is 0
}
 
void draw() {
  if (frameCount % 4 == 0) {
    background(255);
    //lights();
    translate(width/2,2*height/3,-width);
    rotateX(PI/2);
    rotateZ(PI/4);
    for (int i=0; i<N; i++) {
      for (int j=0; j<N; j++) {
        float mapped = map((float)vals[i][j], -.1, .1, 0, 255);
        fill(mapped, 0, 255-mapped, .5*255);
        strokeWeight(1);
        stroke(0);
        //noStroke();
        if (i > 0 & j > 0) {
          beginShape();
          vertex((i-1)*width/N, (j-1)*height/N, (float)vals[i-1][j-1]*height);
          vertex((i)*width/N, (j-1)*height/N, (float)vals[i][j-1]*height);
          vertex((i)*width/N, (j)*height/N, (float)vals[i][j]*height);
          vertex((i-1)*width/N, (j)*height/N, (float)vals[i-1][j]*height);
          endShape(CLOSE);
        }
      }
    }
  }
  
  //step forward in time
  //deepcopy(vals, newvals);
  for (int i=0; i<N; i++) {    
    for (int j=0; j<N; j++) { 
      if (i == 0 | j == 0) {
        newvals[i][j] = 0;
      } else if (i == N-1 | j == N-1) {
        newvals[i][j] = 0;
      } else {
        double ddxx_ddyy = (vals[i+1][j] -2*vals[i][j] + vals[i-1][j])/sq((float)dx)
            + (vals[i][j+1] -2*vals[i][j] + vals[i][j-1])/sq((float)dy);
        double dz = vals[i][j] - oldvals[i][j];
        newvals[i][j] = 2*vals[i][j] - oldvals[i][j] + sq((float)dt)*sq((float)c)*ddxx_ddyy + sq((float)dt)*pressure(t)*dx*dy - b*dz;
      }
    }
  }
  deepcopy(vals, oldvals);
  deepcopy(newvals, vals);

  t = t+dt;
  text((int)t + " seconds", 10, 10);
  text("c = " + c, 10, 25);
  text("dt = " + dt, 10, 40);
  text("P = " +pressure(t), 10, 55);
  
}