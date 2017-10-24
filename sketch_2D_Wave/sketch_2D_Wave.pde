int N = 30;
float t = 0;
float dx = 1/(float)N;
float dy = 1/(float)N;
float dt = .1;
float c = .9;
float b = 0;
float[][] vals = new float[N][N];
float[][] newvals = new float[N][N];
float[][] dvals = new float[N][N];
float[][] newdvals = new float[N][N];

float force(float t) {
  return 10*cos(t/(PI*12));
}
 
void setup() {
  frameRate(100);
  c = c*dx*dy/dt;
  size(500, 500, P3D);
  smooth();
  for (int i=0; i<N; i++) {
    for (int j=0; j<N; j++) {
      vals[i][j] = 0;
      dvals[i][j] = 0;
    }
  }
}
 
void draw() {
  background(255);
  translate(0,2*height/3,-width);
  rotateX(PI/2);
  for (int i=0; i<N; i++) {
    for (int j=0; j<N; j++) {
      float mapped = map(vals[i][j], -height/3, height/3, 0, 255);
      //strokeWeight(5);
      //stroke(mapped, 0, 255-mapped);
      //point(i*width/N,j*height/N, vals[i][j]);
      
      strokeWeight(1);
      fill(mapped, 0, 255-mapped);
      stroke(0);
      if (i > 0) {
        line(
          (i-1)*width/N, (j)*height/N, vals[i-1][j], 
          i*width/N, j*height/N, vals[i][j]
        );
      }
      if (j > 0) {
        line(
          (i)*width/N, (j-1)*height/N, vals[i][j-1], 
          i*width/N, j*height/N, vals[i][j]
        );
      }
      if (i > 0 & j > 0) {
        beginShape();
        vertex((i-1)*width/N, (j-1)*height/N, vals[i-1][j-1]);
        vertex((i)*width/N, (j-1)*height/N, vals[i][j-1]);
        vertex((i)*width/N, (j)*height/N, vals[i][j]);
        vertex((i-1)*width/N, (j)*height/N, vals[i-1][j]);
        endShape(CLOSE);
      }
    }
  }
  
  arrayCopy(vals, newvals);
  arrayCopy(dvals, newdvals);
  for (int i=0; i<N; i++) {    
    for (int j=0; j<N; j++) { 
      if (i == 0 | j == 0) {
        newdvals[i][j] = 0;
        newvals[i][j] = 0;
      } else if (i == N-1 | j == N-1) {
        newdvals[i][j] = 0;
        newvals[i][j] = 0;
      } else {
        float ddvals = sq(c)*(
            (vals[i+1][j] -2*vals[i][j] + vals[i-1][j])/sq(dx)
            + (vals[i][j+1] -2*vals[i][j] + vals[i][j-1])/sq(dy)
          ) - b*dvals[i][j] + force(t)*dx*dy;
        newdvals[i][j] += ddvals*dt;
        newvals[i][j] += newdvals[i][j]*dt;
      }
    }
  }
  arrayCopy(newvals, vals);
  arrayCopy(newdvals, dvals);
  t = t+dt;
}