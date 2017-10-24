int N = 20;
float t = 0;
float dx = 1/(float)N;
float dy = 1/(float)N;
float dt;
float c = 1;
float b = 0;
float[][] vals = new float[N][N];
float[][] newvals = new float[N][N];
float[][] oldvals = new float[N][N];

float force(float t) {
  return 0;// 100*cos(2*PI*t/100);
}
 
void setup() {
  dt = .1*1/c * 1/sqrt(1/sq(dx) + 1/sq(dy));
  println(dt);
  float Cx = c*dt/dx;
  float Cy = c*dt/dy;
  println(sq(Cx)+sq(Cy));
  frameRate(500);
  size(500, 500, P3D);
  smooth();
  for (int i=0; i<N; i++) {
    for (int j=0; j<N; j++) {
      vals[i][j] = pow(2, -(pow((i-N/2),2) + pow((j-N/2),2)));
    }
  }
  arrayCopy(vals, oldvals); //Enforce derivative at t=0 is 0
}
 
void draw() {
  if (frameCount % 4 == 0) {
    background(255);
    translate(0,2*height/3,-width);
    rotateX(PI/2);
    for (int i=0; i<N; i++) {
      for (int j=0; j<N; j++) {
        float mapped = map(vals[i][j], -1, 1, 0, 255);
        fill(mapped, 0, 255-mapped, .7*255);
        strokeWeight(2);
        stroke(0);
        if (i > 0 & j > 0) {
          beginShape();
          vertex((i-1)*width/N, (j-1)*height/N, vals[i-1][j-1]*height);
          vertex((i)*width/N, (j-1)*height/N, vals[i][j-1]*height);
          vertex((i)*width/N, (j)*height/N, vals[i][j]*height);
          vertex((i-1)*width/N, (j)*height/N, vals[i-1][j]*height);
          endShape(CLOSE);
        }
      }
    }
  }
  
  //step forward in time
  for (int i=0; i<N; i++) {    
    for (int j=0; j<N; j++) { 
      if (i == 0 | j == 0) {
        newvals[i][j] = 0;
      } else if (i == N-1 | j == N-1) {
        newvals[i][j] = 0;
      } else {
        float ddxx_ddyy = (vals[i+1][j] -2*vals[i][j] + vals[i-1][j])/sq(dx)
            + (vals[i][j+1] -2*vals[i][j] + vals[i][j-1])/sq(dy);
     
        newvals[i][j] = 2*vals[i][j] - oldvals[i][j] + sq(dt)*sq(c)*ddxx_ddyy + sq(dt)*force(t)*dx*dy;
      }
    }
  }
  arrayCopy(vals, oldvals);
  arrayCopy(newvals, vals);
  t = t+dt;
  text((int)t + " seconds", 10, 10);
  text("c = " + c, 10, 25);
  text("dt = " + dt, 10, 40);
  text("P = " +force(t), 10, 55);
  
}