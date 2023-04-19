using System;
using System.Collections.Generic;
using System.Data.OleDb;
using System.Diagnostics.Eventing.Reader;
using System.Drawing;
using System.Drawing.Imaging;
using System.Linq;
using System.Numerics;
using System.Runtime.InteropServices;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;

namespace ContourDrawing
{
  internal struct Point
  {
    public Point(double x, double y, double z = 0)
    {
      this.x = x;
      this.y = y;
      this.z = z;
    }
    public double x, y, z;
    public bool equals(Point other)
    {
      return (Math.Abs(this.x - other.x) < 0.00000001 && Math.Abs(this.y - other.y) < 0.00000001);
    }
    public double distanceTo(Point other)
    {
      if (!this.equals(other))
      {
        return Math.Pow((Math.Pow(this.x - other.x, 2) + Math.Pow(this.y - other.y, 2)), 0.5);
      }
      else { return 0; }
    }
  }

  internal struct Triangle
  {
    public int p0, p1, p2, t0, t1, t2;
  }

  internal class Program1
  {
    static void Main()
    {
      Application.Run(new Form1());
    }
    public partial class Form1 : Form //I guess this is all necessary?
    {
      public Form1()
      {
        InitializeComponent();
      }
      private void InitializeComponent() //finally
      {
        int height = 767;
        int width = 777;
        double staticContourHeight = 25; //height difference between each contour band on static terrain image
        double staticContourSize = .4; //width of each static contour band
        double dynamicContourSize = 3; //width of the mouse contour band (negative values go down, positive go up)
        double maxZ = 0;
        double[,] depthMap = new double[width, height];
        bool[,] staticContourPixels = new bool[width, height];
        bool[,] scendedPixels = new bool[width, height];
        bool isMoused = false; //well? is it?

        Bitmap backImg = new Bitmap(width, height); //background bitmap
        Bitmap dynamicImg = new Bitmap(width, height);//bitmap containing drawn stuff

        PictureBox displayWindow = new PictureBox(); //container
        displayWindow.Size = new Size(width, height);
        displayWindow.Location = new System.Drawing.Point(0, 0);
        displayWindow.MouseMove += new MouseEventHandler(Moused);
        displayWindow.MouseLeave += new EventHandler(Unmoused);
        displayWindow.Paint += new PaintEventHandler(Painted);
        this.Controls.Add(displayWindow);
        this.Size = new Size(width + 16, height + 39);
        this.CenterToScreen();                               //set up the form

        DrawBmp(backImg); //generate the background image

        




        double[,] testMap = { { 0, 2 }, {0, 2} };
        //Console.Write(testMap[0, 0].ToString() + " " + testMap[1, 0].ToString() + " " 
        //  + testMap[0, 1].ToString() + " " + testMap[1, 1] + "\n") ;

        //var TIN = exportTIN(ref testMap);
        //Triangle[] TINtriangles = TIN.triangles;
        //Point[] TINpoints = TIN.points;
        //List<int> trianglesAtPoint = calcTriangleIndexes(1.7, 0.5, ref testMap);
        Point 
          test1 = new Point(340,50),
          test2 = new Point(200,40), 
          test3 = new Point(300, 500), 
          test4 = new Point(305,0),
          testIntercept = findIntercept(test1,test2,test3,test4);



        Console.WriteLine("test intercept: {0},{1}; {2}", testIntercept.x, testIntercept.y, testIntercept.z);



        (Point[] points, Triangle[] triangles) exportTIN(ref double[,] inputRaster)
        {
          int rasterWidth = inputRaster.GetLength(0);
          int rasterHeight = inputRaster.GetLength(1);
          int numberOfOriginalPoints = (rasterWidth * rasterHeight);
          int numberOfAveragedPoints = (rasterWidth - 1) * (rasterHeight - 1);
          int numberOfPoints = numberOfOriginalPoints + numberOfAveragedPoints;//feel like this is self descriptive.
          int numberOfTriangles = numberOfAveragedPoints * 4;              //bunch of points, originals and renderblock centerpoints
          Point[] points = new Point[numberOfPoints];
          Triangle[] triangles = new Triangle[numberOfTriangles];
          for (int i = 0; i < numberOfPoints; i++)
          {
            points[i] = returnPoint(i, ref inputRaster, numberOfOriginalPoints, numberOfAveragedPoints);//just ask 
          }
          for (int i = 0; i < numberOfTriangles; i++)
          {
            triangles[i] = returnTriangle(i, ref inputRaster); //for the points and triangles
          }
          return (points, triangles);
        }
        Point returnPoint(int pointIndex, ref double[,] inputRaster, int numberOfOriginalPoints = 0, int numberOfAveragedPoints = 0)
        {
          int rasterWidth = inputRaster.GetLength(0);
          int rasterHeight = inputRaster.GetLength(1); //store width and height
          double thisx = 0, thisy = 0, thisz = 0;
          if (numberOfAveragedPoints == 0) //calc numbers if necessary
          {
            numberOfOriginalPoints = (rasterWidth * rasterHeight);
            numberOfAveragedPoints = (rasterWidth - 1) * (rasterHeight - 1);
          }
          if (pointIndex < numberOfOriginalPoints) //if we're in the original point set
          {
            thisx = pointIndex % rasterWidth; //simple enough
            thisy = (int)(pointIndex / rasterWidth);
            thisz = inputRaster[(int)thisx, (int)thisy]; //already known value
          }
          else if (pointIndex < numberOfOriginalPoints + numberOfAveragedPoints)//otherwise we're in the averaged point set
          {
            int lowerx = (int)((pointIndex - numberOfOriginalPoints) % (rasterWidth - 1));//the x and y below the average point
            int lowery = (int)((pointIndex - numberOfOriginalPoints) / (rasterWidth - 1));
            thisx = lowerx + 0.5; //halfway above those
            thisy = lowery + 0.5;
            thisz = (inputRaster[lowerx, lowery] + inputRaster[lowerx + 1, lowery] +
              inputRaster[lowerx, lowery + 1] + inputRaster[lowerx + 1, lowery + 1]) / 4;
          } //z is the average of the four depth values around it
          return new Point { x = thisx, y = thisy, z = thisz };
        }
        Triangle returnTriangle(int triangleIndex, ref double[,] inputRaster)
        {
          int rasterWidth = inputRaster.GetLength(0); //overall width and height
          int rasterHeight = inputRaster.GetLength(1);
          int numberOfOriginalPoints = (rasterWidth * rasterHeight); //number of points in original set
          int numberOfAveragedPoints = (rasterWidth - 1) * (rasterHeight - 1);
          if (triangleIndex > numberOfAveragedPoints * 4 - 1)
          { return new Triangle { p0 = -1, p1 = -1, p2 = -1, t0 = -1, t1 = -1, t2 = -1 }; }//ughhhhh

          int p0 = 0, p1 = 0, p2 = 0, t0 = 0, t1 = 0, t2 = 0; //containers for point and triangle indexes

          int renderBlock = (int)triangleIndex / 4;//block number
          int renderBlockX = renderBlock % (rasterWidth - 1);//renderblock index, and x y
          int renderBlockY = (int)renderBlock / (rasterWidth - 1);

          int blockPoint1 = renderBlock + renderBlockY; //the four points that make up the render block
          int blockPoint2 = blockPoint1 + 1;
          int blockPoint3 = renderBlock + (rasterWidth + 1) + renderBlockY;
          int blockPoint4 = blockPoint3 + 1;

          p1 = renderBlock + numberOfOriginalPoints; //set p1 to averaged center point

          if (triangleIndex % 4 == 0)
          {//first triangle per block
            p0 = blockPoint1;
            p2 = blockPoint2;
            t0 = triangleIndex + 1; //other adjacencies in the block
            t1 = triangleIndex + 3;
            if (renderBlockY > 0) //if not at the top row of render blocks
            {
              t2 = triangleIndex - ((rasterWidth - 1) * 4 - 2); //other triangle is above
            }
            else
            {
              t2 = -1;
            }
          }
          if (triangleIndex % 4 == 1) //second in renderblock
          {
            p0 = blockPoint3;
            p2 = blockPoint1;
            t0 = triangleIndex + 1;
            t1 = triangleIndex - 1;
            if (renderBlockX > 0)
            {
              t2 = triangleIndex - 2;
            }
            else
            {
              t2 = -1;
            }
          }
          if (triangleIndex % 4 == 2) //third
          {
            p0 = blockPoint4;
            p2 = blockPoint3;
            t0 = triangleIndex + 1;
            t1 = triangleIndex - 1;
            if (renderBlockY < (rasterHeight - 2))
            {
              t2 = triangleIndex + ((rasterWidth - 1) * 4 - 2);
            }
            else
            {
              t2 = -1;
            }
          }
          if (triangleIndex % 4 == 3) //fourth and final
          {
            p0 = blockPoint2;
            p2 = blockPoint4;
            t0 = triangleIndex - 3;
            t1 = triangleIndex - 1;
            if (renderBlockX < (rasterWidth - 2))
            {
              t2 = triangleIndex + 2;
            }
            else
            {
              t2 = -1;
            }
          }
          return new Triangle { p0 = p0, p1 = p1, p2 = p2, t0 = t0, t1 = t1, t2 = t2 };
        }
        List<int> calcTriangleIndexes(Point p, ref double[,] inputRaster) //determine which triangles are on a given point
        {
          double x = p.x;
          double y = p.y;
          double distX = x - (int)x;
          double distY = y - (int)y;//how far x&y are from the next lowest whole value
          if (distX > 0.9999999) { x = (int)x + 1; distX = 0; }
          if (distX < 0.00000001) { x = (int)x; distX = 0; }
          if (distY > 0.9999999) { y = (int)y + 1; distY = 0; }//if near an edge, snap to it
          if (distY < 0.00000001) { y = (int)y; distY = 0; }

          List<int> indexes = new List<int> { }; //the indexes of triangles which are at the points

          int rasterWidth = inputRaster.GetLength(0); //overall width and height
          int rasterHeight = inputRaster.GetLength(1);
          if (x > rasterWidth - 1 || x < 0 || y < 0 || y > rasterHeight - 1) { return new List<int> { -1 }; }
          int numberOfOriginalPoints = (rasterWidth * rasterHeight); //number of points in original set
          int numberOfAveragedPoints = (rasterWidth - 1) * (rasterHeight - 1);
          int firstTriInBlock = ((int)y * (rasterWidth - 1) + (int)x) * 4;//index of lowest number triangle in render block

          if (distX == 0 && distY == 0)
          {                             //we're at a whole number point
            int triURR = firstTriInBlock - ((rasterWidth - 1) * 4 - 2), triUUR = triURR - 1, triUUL = triUUR - 2, triULL = triUUL - 1,
              triDRR = firstTriInBlock, triDDR = triDRR + 1, triDDL = triDRR - 1, triDLL = triDDL - 3;
            List<int> trisUR = new List<int> { triURR, triUUR };//the 8 triangles around, and blocks to store them
            List<int> trisUL = new List<int> { triUUL, triULL };
            List<int> trisDL = new List<int> { triDLL, triDDL };//worth mentioning that (0,0) is top left of raster for these UDLR
            List<int> trisDR = new List<int> { triDDR, triDRR };
                                                              //check for edges, add triangles if possible
            if (x != rasterWidth - 1 && y > 0) { indexes.AddRange(trisUR); } //can add up right two, if not at right edge or top
            if (x > 0 && y > 0) { indexes.AddRange(trisUL); }
            if (x > 0 && y != rasterHeight - 1) { indexes.AddRange(trisDL); }
            if (x != rasterWidth - 1 && y != rasterHeight - 1) { indexes.AddRange(trisDR); }
            return indexes;
          }

          else if (Math.Abs(distX - 0.5) < 0.00000001 && Math.Abs(distY - 0.5) < 0.00000001)
          {   //we're at an averaged point, add whole renderblock of triangles
            indexes.Add(firstTriInBlock); indexes.Add(firstTriInBlock + 1);
            indexes.Add(firstTriInBlock + 2); indexes.Add(firstTriInBlock + 3);
            return indexes;
          }

          else if (distX < 0.00000001  && distY > 0.00000001) //on a vertical line
          {
            if (x > 0) { indexes.Add(firstTriInBlock - 1); }
            if (x < rasterWidth - 1) { indexes.Add(firstTriInBlock + 1); }
            return indexes;
          }
          else if (distY < 0.00000001  && distX > 0.00000001) //on a horizontal line
          {
            if (y > 0) { indexes.Add(firstTriInBlock - ((rasterWidth - 1) * 4 - 2)); }
            if (y < rasterHeight - 1) { indexes.Add(firstTriInBlock); }
            return indexes;
          }

          else if (distY > 0 && distX > 0)//somewhere else in the render block
          {
            if (Math.Abs(distY - distX) < 0.00000001)//diagonal line UL to DR 
            {
              if (distX < 0.5) { indexes.Add(firstTriInBlock); indexes.Add(firstTriInBlock + 1); }
              else { indexes.Add(firstTriInBlock + 2); indexes.Add(firstTriInBlock + 3); }
            }
            else if (Math.Abs(distY + distX - 1) < 0.00000001)//diagonal line DL to UR
            {
              if (distX < distY) { indexes.Add(firstTriInBlock + 1); indexes.Add(firstTriInBlock + 2); }
              else { indexes.Add(firstTriInBlock); indexes.Add(firstTriInBlock + 3); }
            }
            else //we're not in kansas anymore
            {
              if (distX + distY < 1) //UL half
              {
                if (distX > distY) { indexes.Add(firstTriInBlock); }//U quarter
                else { indexes.Add(firstTriInBlock + 1); }//L quarter
              }
              else  //DR half
              {
                if (distX > distY) { indexes.Add(firstTriInBlock + 3); }//R quarter
                else { indexes.Add(firstTriInBlock + 2); }//D quarter
              }
            }
            return indexes;
          }
          return indexes;//again, with feeling
        }




        void Moused(object sender, MouseEventArgs e) //lol get moused
        {
          isMoused = true; //sure is

          Point mouseP = new Point();
          int mouseX = e.X;
          int mouseY = e.Y;
          double mouseHeight = depthMap[mouseX, mouseY];
          mouseP.x = mouseX; mouseP.y = mouseY; mouseP.z = mouseHeight;

          Console.WriteLine(mouseHeight);//just for funsies

          drawContour(ref scendedPixels, mouseHeight, dynamicContourSize);//contour at mouseheight

          Parallel.For(0, 2, i =>
          { //do these in seperate threads so it's faster
            if (i == 0)
            {
              Scender(mouseP, 1, ref depthMap); //call Ascend
            }
            else
            {
              Scender(mouseP, -1, ref depthMap);//call Descend
            }
          });
          updateScreen(mouseHeight); //generate the dynamic bitmap and paint it to the display
        }
        void Scender(Point currentPoint, int up, ref double[,] inputRaster) //Ascend or Descend recursively with all new magic triangles
        {
          double x = currentPoint.x;
          double y = currentPoint.y;
          scendedPixels[ (int)x, (int)y ] = true; //color current pixel

          if (x <= 0 || x > inputRaster.GetLength(0) || y <= 0 || y > inputRaster.GetLength(1)) { return; } //end if edge
          Point steepestPoint = currentPoint;

          List<int> tris = calcTriangleIndexes(currentPoint, ref inputRaster); //for each available triangle at a point
          //Console.WriteLine(tris.Count());
          foreach (int triIndex in tris)
          {
            Point triBestPoint = findSteepestPoint(triIndex, currentPoint, up, ref inputRaster);//find steepest point
            double triBestSlope = findZSlope(currentPoint, triBestPoint);//calculate slope
            double steepestSlope = findZSlope(currentPoint, steepestPoint);//slope of steepest point so far
            
            if (triBestSlope * up > steepestSlope * up) { steepestPoint = triBestPoint; }
            
          }

          if (!steepestPoint.equals(currentPoint))// && up ==1)
          {
            Scender(steepestPoint, up, ref inputRaster);//recur
          } 
          else if (steepestPoint.equals(currentPoint))
          {
            //we made it?
            return;
          }

          //some of this stuff might become relevant for flat areas?
          /*
          int i = 8; //trust the math
          double originDepth = depthMap[X, Y]; //depth of circle origin
          //double originDepth = depth on line between closest pixel center and next pixel over in correct direction 
          bool scended = false; //are we there yet?
          while (scended == false) //no, not yet
          {
            double steepest = originDepth; //holder for greatest inspected value
            int steepestCount = 0; //how many pixels are at the steepest height
            double radius = i / (Math.PI * 2); //backward calculate radius from required circumference. Trust the math
            double[][] inspPixels = new double[i][]; //container for current circle data
            for (int k = 0; k < (i); k++) //ITERATE
            {
              double angle = Math.PI * 2 / i * k; //fraction of tau to check for a pixel
              int inspX = (int)X + (int)(Math.Sin(angle) * radius);//trust
              int inspY = (int)Y + (int)(Math.Cos(angle) * radius);//the Math
              if (height > inspY && inspY >= 0 && width > inspX && inspX >= 0) //in bounds?
              {
                double inspDepth = depthMap[inspX, inspY]; //depth of current pixel
                inspPixels[k] = new double[] { inspX, inspY, inspDepth }; //write current pixel info to array
                if ((inspDepth - originDepth) * up > 0)//diffbetween current & origin depth in right direction
                {
                  scended = true;//WE'RE THERE
                  if ((inspDepth - steepest) * up > 0)//diffbetween current and steepest is larger in right direction
                  {
                    steepest = inspDepth; //record new highest or lowest point
                  }
                }
              }
              else //we're out of bounds so
              {
                inspPixels[k] = new double[] { 0, 0, -1 }; //set the inspected pixel to a negated value
              }
            } //INITIAL ITERATION COMPLETE

            for (int m = 0; m < inspPixels.GetLength(0); m++) //loop through the circle to
            {
              if (inspPixels[m][2] == steepest) //check for pixels at steepest height, and
              {
                steepestCount++; //add up the total number of pixels at that height
              }
            } //steepestcount loop complete, now

            bool PeakOrValley = true;//we wouldn't want to keep climbing at the peak or floor
            int currentSteepest = 0; //how far through the list of steepest are we
            for (int u = 0; u < inspPixels.GetLength(0); u++) //loop through current circle pixels, decide actions
            {
              double[] current = inspPixels[u]; //current pixel
              if (current[2] != -1) //in bounds?
              {
                if (current[2] == steepest && steepest != originDepth)//if pixel is among the steepest
                {
                  currentSteepest++;//this is how far we are through steepest
                  if (scendedPixels[(int)current[0], (int)current[1]] == false)//if pixel not already colored
                  {
                    int doubleSteepest = currentSteepest * 2;
                    int fractionOfSteep = steepestCount - doubleSteepest; //this stuff COULD be used to choose middle path
                    //Console.WriteLine(fractionOfSteep);
                    if (true)//-1 <= fractionOfSteep && fractionOfSteep <= 1)//are we at a middle steepest values
                    {
                      Scender((int)current[0], (int)current[1], up);//recur at pixel
                    }
                  }
                } //whether or not at the steepest
                if (((current[2] - originDepth) * up) >= 0)//current depth level or in right direction
                {
                  scendedPixels[(int)current[0], (int)current[1]] = true;//color pixel at draw time
                  PeakOrValley = false; //can't be
                }
              }
            }
            if (PeakOrValley == true || originDepth == 0)
            {
              scended = true;
            }
            i += 4; //trust the math
          }*/
        }
        Point findSteepestPoint(int triIndex, Point startPoint, int up, ref double[,] inputRaster)
        {

          int rasterWidth = inputRaster.GetLength(0); //overall width and height
          int rasterHeight = inputRaster.GetLength(1);
          int numberOfOriginalPoints = (rasterWidth * rasterHeight); //number of points in original set
          int numberOfAveragedPoints = (rasterWidth - 1) * (rasterHeight - 1);
          if (triIndex >= 0 && triIndex < numberOfOriginalPoints + numberOfAveragedPoints)
          {
            Triangle derefTri = returnTriangle(triIndex, ref inputRaster);
            Point pA = returnPoint(derefTri.p1, ref inputRaster); //always middle point
            Point pB = returnPoint(derefTri.p0, ref inputRaster);
            Point pC = returnPoint(derefTri.p2, ref inputRaster);

            Point normalVector = findNormal(pA, pB, pC);//<take the normal \/make it relative to startpoint
            Point normalVectorPoint = new Point { x = normalVector.x + startPoint.x, y = normalVector.y + startPoint.y, z = normalVector.z + startPoint.z };
            if (normalVectorPoint.equals(startPoint))
            {                                         //opposite normal if same as startpoint
              normalVector = findNormal(pA, pC, pB);
              normalVectorPoint = new Point { x = normalVector.x + startPoint.x, y = normalVector.y + startPoint.y, z = normalVector.z + startPoint.z };
            }

            Point[] intercepts = { //fill an array with the intersections of the steepest(normal) line and the triangle lines
            findIntercept(startPoint, normalVectorPoint, pA, pB),
            findIntercept(startPoint, normalVectorPoint, pB, pC),//if outside bounds of triangle, z==-1
            findIntercept(startPoint, normalVectorPoint, pC, pA)
            };
            //Console.WriteLine("\nStartPoint: {0},{1}; {2}\n", startPoint.x, startPoint.y, startPoint.z);
            foreach (Point intercept in intercepts)
            {
              //Console.Write("intercept: {0},{1};{2}, slope {3}\n", intercept.x, intercept.y, intercept.z, findZSlope(startPoint, intercept));
              if (intercept.equals(startPoint))
              {
                //startPoint.z = intercept.z;
              }
            }
            Point triBestPoint = startPoint;

            foreach (Point intercept in intercepts)
            {
              double intSlope = findZSlope(startPoint, intercept); // slope from start to inspected intercept
              double bestSlope = 0; 
              if (!triBestPoint.equals(startPoint)) { bestSlope = findZSlope(startPoint, triBestPoint); } 

              if (intercept.z != -1 && intSlope * up > bestSlope * up)//steeper in correct direction
              {
                triBestPoint = intercept;
              }
            }
            if (triBestPoint.equals(startPoint))//if no good intercepts
            {
              double slopeA = findZSlope(startPoint, pA); //examine the three points
              double slopeB = findZSlope(startPoint, pB);
              double slopeC = findZSlope(startPoint, pC);
              if (slopeA * up > 0) { triBestPoint = pA; }
              if (slopeB * up > 0 && slopeB * up > slopeA * up) { triBestPoint = pB; }//pick the best one
              if (slopeC * up > 0 && slopeC * up > slopeB * up && slopeC * up > slopeA * up) { triBestPoint = pC; }
            }
            return triBestPoint;
          }
          else { return startPoint; }
        }
        Point findIntercept(Point L1P1, Point L1P2, Point L2P1, Point L2P2)
        {
          double x = -1, y = -1, z = -1;//default/error values

          double L1Rise = L1P1.y - L1P2.y;//rise?
          double L1Run = L1P1.x - L1P2.x;//run?
          double L1M; if (L1Run != 0) { L1M = L1Rise / L1Run; } else { L1M = double.PositiveInfinity; }
          double L1B = L1P1.y - L1M * L1P1.x;

          double L2Rise = L2P1.y - L2P2.y;
          double L2Run = L2P1.x - L2P2.x;
          double L2M; if (L2Run != 0) { L2M = L2Rise / L2Run; } else { L2M = double.PositiveInfinity; }
          double L2B = L2P1.y - L2M * L2P1.x;

          if (L1Run == 0 && L2Run != 0) //first line vertical
          {
            x = L1P1.x;
            y = L2M * x + L2B;
          }
          else if (L2Run == 0 && L1Run != 0) //second line vertical
          {
            x = L2P1.x;
            y = L1M * x + L1B;
          }
          else // calc intercept
          {
            if (L1Run != L2Run) //not parallel
            {
              x = (L2B - L1B)/(L1M - L2M);
              y = L1M * x + L1B;
            } else { return new Point(x, y, z); }
          }
          if (L2P1.x >= x && x >= L2P2.x || L2P1.x <= x && x <= L2P2.x 
           && L2P1.y >= y && y >= L2P2.y || L2P1.y <= y && y <= L2P2.y
           && L1P1.y == L1P1.x * L2M + L2B)//intercept is on and within L2 segment
          {
            Point intercept = new Point { x = x, y = y, z = -1 };
            double distFromP1 = L2P1.distanceTo(intercept);
            double distFromP2 = L2P2.distanceTo(intercept);
            double L2Length = distFromP1 + distFromP2;
            double fractFromP1 = distFromP1 / L2Length;
            double zDiff = L2P2.z - L2P1.z;
            z = zDiff * fractFromP1 + L2P1.z;
          }
          return new Point { x = x, y = y, z = z };
        }
        double findHeight(Point p, ref double[,] inputRaster) //return height of point in TIN
        {
          List<int> trisAtPoint = calcTriangleIndexes(p, ref inputRaster);
          Triangle tri = returnTriangle(trisAtPoint[0], ref inputRaster);
          Point pA = returnPoint(tri.p1, ref inputRaster); 
          Point pB = returnPoint(tri.p0, ref inputRaster);
          Point pC = returnPoint(tri.p2, ref inputRaster);
          Point normal = findNormal(pA, pB, pC);
          double A = normal.x, B = normal.y, C = normal.z, D = -(A * pA.x) - (B * pA.y) - (C * pA.z); //plane equation
          return -(A * p.x + B * p.y + D) / C;
        }
        Point findNormal(Point pA, Point pB, Point pC)//do it
        {

          Point vectAB = new Point { x = pB.x - pA.x, y = pB.y - pA.y, z = pB.z - pA.z };
          Point vectAC = new Point { x = pC.x - pA.x, y = pC.y - pA.y, z = pC.z - pA.z };
          return new Point
          {
            x = vectAB.y * vectAC.z - vectAB.z * vectAC.y,
            y = vectAB.z * vectAC.x - vectAB.x * vectAC.z,
            z = vectAB.x * vectAC.y - vectAB.y * vectAC.x
          };
        }
        double findZSlope(Point p1, Point p2)//ezpz
        {
          double zDiff = p2.z - p1.z;
          double dist = p1.distanceTo(p2);
          if (dist > 0.00000000000001) { return zDiff / dist; }
          else { return 0; }
        }



        void updateScreen(double mouseHeight = 0) //mouse event draw, creates an overlay bmp with contour & up/down line
        {
          var byteMe = bmp2Array(dynamicImg);//do it
          byte[] byteData = byteMe.byteData;
          BitmapData rawData = byteMe.rawData;//bit silly, innit
          Parallel.For(0, (rawData.Stride * rawData.Height / 4), i => //loop through the rows of data, in parallel
          {
            int iter = i * 4; //4 bytes per pixel
            int col = (iter % rawData.Stride) / 4;//calc pixel x value
            int row = iter / rawData.Stride; //calc pixel y value
            if (scendedPixels[col, row] == true)
            {       //if pixel in scender map, color the pixel
              byteData[iter] = (byte)(depthMap[col, row]);        //blue
              byteData[iter + 1] = (byte)(depthMap[col, row] * .4);    //green
              byteData[iter + 2] = (byte)(-depthMap[col, row] + 255);  //red
              byteData[iter + 3] = 200;                                //alpha
            }
            else //otherwise,
            {
              byteData[iter + 3] = 0; //make it transparent
            }
          });
          scendedPixels = new bool[width, height]; //clear the scended array
          bmpFromArray(byteData, rawData, dynamicImg); //return data to bitmap, from the array
          displayWindow.Refresh(); //paint!
        }
        (byte[] byteData, BitmapData rawData) bmp2Array(Bitmap bmp) //convert bitmap into byte array & bitmapdata
        {
          BitmapData rawData = bmp.LockBits(new Rectangle(0, 0, width, height), ImageLockMode.ReadWrite, bmp.PixelFormat);
          int dataSize = rawData.Stride * rawData.Height; //Actually not sure
          byte[] byteData = new byte[dataSize];           //if I'm being honest,
          System.Runtime.InteropServices.Marshal.Copy(rawData.Scan0, byteData, 0, dataSize);//I basically stole this
          return (byteData, rawData); //pretty cool you can return multiple values
        }
        void bmpFromArray(byte[] byteData, BitmapData rawData, Bitmap bmp) //copy data from byte array back to bmp
        {
          System.Runtime.InteropServices.Marshal.Copy(byteData, 0, rawData.Scan0, byteData.Length);//stole this too
          bmp.UnlockBits(rawData); //BE FREEEEE
        }

        void Unmoused(object sender, EventArgs e) //no mouse 4 u
        {
          isMoused = false; //not anymore
          displayWindow.Refresh();
        }

        void Painted(object sender, PaintEventArgs p) //Ponting
        {
          Graphics g = p.Graphics;
          g.DrawImageUnscaled(backImg, 0, 0); //something to do with painting
          if (isMoused == true) //is it?
          {
            g.DrawImageUnscaled(dynamicImg, 0, 0);//mouse that shit
          }
        }

        void DrawBmp(Bitmap bmp) //initializing draw, creates background image of contoured terrain
        {
          var byteMe = (bmp2Array(bmp)); // >:3
          byte[] byteData = byteMe.byteData;
          BitmapData rawData = byteMe.rawData;
          Parallel.For(0, (rawData.Stride * rawData.Height / 4), i => //parallel loop through the pixel rows, generate depths
          {
            int iter = i * 4;
            int col = (iter % rawData.Stride) / 4; //calc x value
            int row = iter / rawData.Stride; //calc y value
            double colorval = DepthFinder(col, row); //calculate depth of pixel
            depthMap[col, row] = colorval; //save depth in array
            if (colorval > maxZ)
            {
              maxZ = colorval;
            }
          });

          double scaleFactor = 255 / maxZ;
          for (int i = 0; i < depthMap.GetLength(0); i++)
          {
            for (int j = 0; j < depthMap.GetLength(1); j++)
            {
              depthMap[i, j] = scaleFactor * depthMap[i, j];
            }
          }

          Parallel.For(1, (int)(maxZ / staticContourHeight + 1), i =>
          {
            drawContour(ref staticContourPixels, staticContourHeight * i, staticContourSize);
          });

          Parallel.For(0, (rawData.Stride * rawData.Height / 4), i =>
          {
            int iter = i * 4;
            int col = (iter % rawData.Stride) / 4; //calc x value
            int row = iter / rawData.Stride; //calc y value
            byteData[iter + 3] = 255; //set alpha to full
            for (int k = iter; k < iter + 3; k++) //set color values 
            {
              if (staticContourPixels[col, row])//colorval % staticContourHeight == 0)
              {
                byteData[k] = 0; //set black for contour
              }
              else
              {
                byteData[k] = (byte)depthMap[col, row]; //set grayscale as height
              }
            }
          });
          bmpFromArray(byteData, rawData, bmp);
        }
        double DepthFinder(int x, int y) //calculates & returns the depth of a given pixel, from its x,y coordinate
        {
          
          return (double)(1 * y + 1 * x) / 5;

          double mxRow = ((double)(y + 1) / height); //fraction of full width or height
          double mxCol = ((double)(x + 1) / width);
          double tauRow = mxRow * Math.PI * 2; //fraction times tau
          double tauCol = mxCol * Math.PI * 2;
          double sinRow = Math.Sin(tauRow); //sin of tau at that angle
          double sinCol = Math.Sin(tauCol);
          double somnfancy = 120 * Math.Abs((((sinRow) * (sinCol)) + 1) / 2 - (sinRow - sinCol)); //WOOOOOOOO
          //return somnfancy;

        }


        void drawContour(ref bool[,] binaryArray, double contourHeight, double bandSize)
        {
          bool[,] localArray = binaryArray;//local container for ref array so lambda can use it
          Parallel.For(0, height, y =>
          {
            for (int x = 0; x < width; x++) //loop through pixels
            {
              double currentHeight = depthMap[x, y];//height of center pixel

              if ((currentHeight <= contourHeight) && (currentHeight >= contourHeight + bandSize) && (bandSize < 0) 
                  ||
                  (currentHeight >= contourHeight) && (currentHeight <= contourHeight + bandSize) && (bandSize > 0))
              {
                localArray[x, y] = true; //if pixel in the band, color it
              }

              if ((y + 1) < height) //if not out of bounds, check for contours between adjacent pixels
              {
                colorClosestPixel(ref localArray, x, y, x, (y + 1), contourHeight, bandSize);
              }
              if ((x + 1) < width)
              {
                colorClosestPixel(ref localArray, x, y, (x + 1), y, contourHeight, bandSize);
              }

            }
          });

        }
        void colorClosestPixel(ref bool[,] binaryArray, int x1, int y1, int x2, int y2, double contourHeight, double bandSize)
        {
          double currentHeight = depthMap[x1, y1]; //heights of two pixels to check
          double inspHeight = depthMap[x2, y2];
          if (bandSize < Math.Abs(currentHeight - inspHeight)) //if the band is smaller than the gap
          {
            //if the two contour heights are between the inspected heights
            if (((inspHeight >= contourHeight && currentHeight <= contourHeight) ||
                (inspHeight <= contourHeight && currentHeight >= contourHeight)) 
                &&
                ((inspHeight >= contourHeight - bandSize && currentHeight <= contourHeight - bandSize) ||
                (inspHeight <= contourHeight - bandSize && currentHeight >= contourHeight - bandSize)))
            {
              //if the closest distance to insp is greater than the closest to current
              if (Math.Min(Math.Abs(inspHeight - contourHeight), Math.Abs(inspHeight - (contourHeight - bandSize))) >
                  Math.Min(Math.Abs(currentHeight - contourHeight), Math.Abs(currentHeight - (contourHeight - bandSize))))
              {
                binaryArray[x1, y1] = true; //color current
              }
              else //otherwise the closest is insp
              {
                binaryArray[x2, y2] = true; //color insp
              }
            }
          }
        }
      }
    }
  }
}