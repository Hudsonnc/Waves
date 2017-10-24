int N = 100;
float t = 0;
float dx = 1/(float)N;
float dt = .1;
float c = .9;
float b = 0;
float[] vals = new float[N];
float[] newvals = new float[N];
float[] dvals = new float[N];
float[] newdvals = new float[N];

float force(float t) {
  return 1000*cos(t/(PI/2));
}
 
void setup() {
  frameRate(100);
  c = c*dx/dt;
  size(500, 500, P3D);
  strokeWeight(20);
  smooth();
  for (int i=0; i<vals.length; i++) {
    vals[i] = height/2;
    dvals[i] = 0;
  }
}
 
void draw() {
  background(255);
  //rectMode(CENTER);
  translate(0,700,-600);
  rotateX(PI/2);
  for (int i=0; i<N; i++) {
    for (int j=0; j<N; j++) {
      float mapped = map(vals[i], height/3, 2*height/3, 0, 255);
      
      stroke(mapped, 0, 255-mapped);
      point(i*width/N,j*height/N, vals[i]);
    }
  }
  //print(frameRate);
  
  
  arrayCopy(vals, newvals);
  arrayCopy(dvals, newdvals);
  for (int i=0; i<vals.length; i++) {    
    if (i == 0) {
      newdvals[i] = 0;
      newvals[i] = height/2;
    } else if (i == vals.length-1) {
      newdvals[i] = 0;
      newvals[i] = height/2;
    } else {
      float ddvals = sq(c)*(vals[i+1] -2*vals[i] + vals[i-1])/sq(dx) - b*dvals[i] + force(t)*dx;
      newdvals[i] += ddvals*dt;
      newvals[i] += newdvals[i]*dt;
    }
  }
  arrayCopy(newvals, vals);
  arrayCopy(newdvals, dvals);
  t = t+dt;
}