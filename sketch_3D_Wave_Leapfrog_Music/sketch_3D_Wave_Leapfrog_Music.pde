import ddf.minim.*;
 
Minim minim;
AudioPlayer song;

int N = 20;
double t = 0;
double dx = 1/(double)N;
double dy = 1/(double)N;
double dz = 1/(double)N;
double dt = 1/240.0;
double c = .1;
double b = .01;
double[][][] vals = new double[N][N][N];
double[][][] newvals = new double[N][N][N];
double[][][] oldvals = new double[N][N][N];

double pressure(double t) {
  return (song.left.get(0) + song.right.get(0))*90000; 
}

void deepcopy(double[][][] src, double[][][] dest) {
  for (int i = 0; i < src.length; i++) {
    for (int j = 0; j < src.length; j++) {
      arrayCopy(src[i][j], dest[i][j]);
    }
  }
}
 
void setup() {
  c = c*.9*1/dt * 1/sqrt(1/sq((float)dx) + 1/sq((float)dy) + 1/sq((float)dz));
  frameRate(240);
  size(1000,1000, P3D);
  smooth();
  for (int i=0; i<N; i++) {
    for (int j=0; j<N; j++) {
      for (int k=0; k<N; k++) {
        vals[i][j][k] = 0; //.5*exp(-.02*(pow(2*(i-N/2),2) + pow((j-N/2),2) + pow((k-N/2),2)));
      }
    }
  }
  deepcopy(vals, oldvals); //Enforce derivative at t=0 is 0
  
  minim = new Minim(this);
 
  // this loads mysong.wav from the data folder
  song = minim.loadFile("talktome.mp3");
  song.play();
  
}
 
void draw() {
  if (frameCount % 2 == 0) {
    background(0);
    translate(width/2,3*height/4,-width);
    rotateX(PI/2);
    rotateZ(PI/4);
    for (int i=0; i<N; i++) {
      for (int j=0; j<N; j++) {
        for (int k=0; k<N; k++) {
          float mapped = map((float)vals[i][j][k], -.1, .1, 0, 255);
          stroke(mapped, 0, 255-mapped, .5*255);
          strokeWeight(mapped/255*10);
          point((i*width)/N*2/3, (j*height)/N*2/3, (k*height)/N*2/3);
        }
      }
    }
  }
  
  //step forward in time
  for (int i=0; i<N; i++) {    
    for (int j=0; j<N; j++) { 
      for (int k=0; k<N; k++) {
        if (i == 0 | j == 0 | k == 0) {
          newvals[i][j][k] = 0;
        } else if (i == N-1 | j == N-1 | k == N-1) {
          newvals[i][j][k] = 0;
        } else {
          double ddxx_ddyy_ddzz = (vals[i+1][j][k] -2*vals[i][j][k] + vals[i-1][j][k])/sq((float)dx)
              + (vals[i][j+1][k] -2*vals[i][j][k] + vals[i][j-1][k])/sq((float)dy)
              + (vals[i][j][k+1] -2*vals[i][j][k] + vals[i][j][k-1])/sq((float)dz);
          double dvals = vals[i][j][k] - oldvals[i][j][k];
          newvals[i][j][k] = 2*vals[i][j][k] - oldvals[i][j][k] + sq((float)dt)*sq((float)c)*ddxx_ddyy_ddzz + sq((float)dt)*pressure(t)*dx*dy*dz - b*dvals;
        }
      }
    }
  }
  deepcopy(vals, oldvals);
  deepcopy(newvals, vals);

  t = t+dt;
  fill(0);
  text((int)t + " seconds", 10, 10);
  text("c = " + c, 10, 25);
  text("dt = " + dt, 10, 40);
  text("P = " +pressure(t), 10, 55);
  
}