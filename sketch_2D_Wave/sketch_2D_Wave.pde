int N = 30;
float t = 0;
float dx = 1/(float)N;
float dy = 1/(float)N;
float dt;
float c = 1;
float b = 1;
float[][] vals = new float[N][N];
float[][] newvals = new float[N][N];
float[][] dvals = new float[N][N];
float[][] newdvals = new float[N][N];

float force(float t) {
  return cos(2*PI*t);
  //return(0);
}
 
void setup() {
  dt = .1* 1/c * 1/sqrt(1/sq(dx) + 1/sq(dy));
  print(dt);
  frameRate(500);
  size(500, 500, P3D);
  smooth();
  for (int i=0; i<N; i++) {
    for (int j=0; j<N; j++) {
      vals[i][j] = 0;//.0000001*pow(2, -pow((i-N/2),2)-pow((j-N/2),2));
      dvals[i][j] = 0;
    }
  }
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
        newvals[i][j] += dvals[i][j]*dt;
      }
    }
  }
  arrayCopy(newvals, vals);
  arrayCopy(newdvals, dvals);
  t = t+dt;
  text((int)t + " seconds", 10, 10);
  text("c = " + c, 10, 25);
  text("dt = " + dt, 10, 40);
  text("P = " +force(t), 10, 55);
  
}