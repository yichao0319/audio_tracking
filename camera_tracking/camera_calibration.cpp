#include "opencv2/core.hpp"
#include <opencv2/core/utility.hpp>
#include "opencv2/imgproc.hpp"
#include "opencv2/calib3d.hpp"
#include "opencv2/imgcodecs.hpp"
#include "opencv2/videoio.hpp"
#include "opencv2/highgui.hpp"

#include <cctype>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <iostream>

using namespace cv;
using namespace std;

const char * usage =
" \nexample command line for calibration from a live feed.\n"
"   calibration  -w 4 -h 5 -s 0.025 -o camera.yml -op -oe\n"
" \n"
" example command line for calibration from a list of stored images:\n"
"   imagelist_creator image_list.xml *.png\n"
"   calibration -w 4 -h 5 -s 0.025 -o camera.yml -op -oe image_list.xml\n"
" where image_list.xml is the standard OpenCV XML/YAML\n"
" use imagelist_creator to create the xml or yaml list\n"
" file consisting of the list of strings, e.g.:\n"
" \n"
"<?xml version=\"1.0\"?>\n"
"<opencv_storage>\n"
"<images>\n"
"view000.png\n"
"view001.png\n"
"<!-- view002.png -->\n"
"view003.png\n"
"view010.png\n"
"one_extra_view.jpg\n"
"</images>\n"
"</opencv_storage>\n";




const char* liveCaptureHelp =
    "When the live video from camera is used as input, the following hot-keys may be used:\n"
        "  <ESC>, 'q' - quit the program\n"
        "  'g' - start capturing images\n"
        "  'u' - switch undistortion on/off\n";

static void help()
{
    printf( "This is a camera calibration sample.\n"
        "Usage: calibration\n"
        "     -w <board_width>         # the number of inner corners per one of board dimension\n"
        "     -h <board_height>        # the number of inner corners per another board dimension\n"
        "     [-pt <pattern>]          # the type of pattern: chessboard or circles' grid\n"
        "     [-n <number_of_frames>]  # the number of frames to use for calibration\n"
        "                              # (if not specified, it will be set to the number\n"
        "                              #  of board views actually available)\n"
        "     [-d <delay>]             # a minimum delay in ms between subsequent attempts to capture a next view\n"
        "                              # (used only for video capturing)\n"
        "     [-s <squareSize>]       # square size in some user-defined units (1 by default)\n"
        "     [-o <out_camera_params>] # the output filename for intrinsic [and extrinsic] parameters\n"
        "     [-op]                    # write detected feature points\n"
        "     [-oe]                    # write extrinsic parameters\n"
        "     [-zt]                    # assume zero tangential distortion\n"
        "     [-a <aspectRatio>]       # fix aspect ratio (fx/fy)\n"
        "     [-p]                     # fix the principal point at the center\n"
        "     [-v]                     # flip the captured images around the horizontal axis\n"
        "     [-V]                     # use a video file, and not an image list, uses\n"
        "                              # [input_data] string for the video file name\n"
        "     [-su]                    # show undistorted images after calibration\n"
        "     [input_data]             # input data, one of the following:\n"
        "                              #  - text file with a list of the images of the board\n"
        "                              #    the text file can be generated with imagelist_creator\n"
        "                              #  - name of video file with a video of the board\n"
        "                              # if input_data not specified, a live view from the camera is used\n"
        "\n" );
    printf("\n%s",usage);
    printf( "\n%s", liveCaptureHelp );
}

enum { DETECTION = 0, CAPTURING = 1, CALIBRATED = 2 };
enum Pattern { CHESSBOARD, CIRCLES_GRID, ASYMMETRIC_CIRCLES_GRID };

static double computeReprojectionErrors(
        const vector<vector<Point3f> >& objectPoints,
        const vector<vector<Point2f> >& imagePoints,
        const vector<Mat>& rvecs, const vector<Mat>& tvecs,
        const Mat& cameraMatrix, const Mat& distCoeffs,
        vector<float>& perViewErrors )
{
    vector<Point2f> imagePoints2;
    int i, totalPoints = 0;
    double totalErr = 0, err;
    perViewErrors.resize(objectPoints.size());

    for( i = 0; i < (int)objectPoints.size(); i++ )
    {
        projectPoints(Mat(objectPoints[i]), rvecs[i], tvecs[i],
                      cameraMatrix, distCoeffs, imagePoints2);
        err = norm(Mat(imagePoints[i]), Mat(imagePoints2), NORM_L2);
        int n = (int)objectPoints[i].size();
        perViewErrors[i] = (float)std::sqrt(err*err/n);
        totalErr += err*err;
        totalPoints += n;
    }

    return std::sqrt(totalErr/totalPoints);
}

static void calcChessboardCorners(Size boardSize, float squareSize, vector<Point3f>& corners, Pattern patternType = CHESSBOARD)
{
    corners.resize(0);

    switch(patternType)
    {
      case CHESSBOARD:
      case CIRCLES_GRID:
        for( int i = 0; i < boardSize.height; i++ )
            for( int j = 0; j < boardSize.width; j++ )
                corners.push_back(Point3f(float(j*squareSize),
                                          float(i*squareSize), 0));
        break;

      case ASYMMETRIC_CIRCLES_GRID:
        for( int i = 0; i < boardSize.height; i++ )
            for( int j = 0; j < boardSize.width; j++ )
                corners.push_back(Point3f(float((2*j + i % 2)*squareSize),
                                          float(i*squareSize), 0));
        break;

      default:
        CV_Error(Error::StsBadArg, "Unknown pattern type\n");
    }
}

static bool runCalibration( vector<vector<Point2f> > imagePoints,
                    Size imageSize, Size boardSize, Pattern patternType,
                    float squareSize, float aspectRatio,
                    int flags, Mat& cameraMatrix, Mat& distCoeffs,
                    vector<Mat>& rvecs, vector<Mat>& tvecs,
                    vector<float>& reprojErrs,
                    double& totalAvgErr)
{
    cameraMatrix = Mat::eye(3, 3, CV_64F);
    if( flags & CALIB_FIX_ASPECT_RATIO )
        cameraMatrix.at<double>(0,0) = aspectRatio;

    distCoeffs = Mat::zeros(8, 1, CV_64F);

    vector<vector<Point3f> > objectPoints(1);
    calcChessboardCorners(boardSize, squareSize, objectPoints[0], patternType);

    objectPoints.resize(imagePoints.size(),objectPoints[0]);


    // ycc
    // vector<vector<Point2f> > imagePoints;
    // int cnt1 = 0;
    // for (vector<vector<Point2f> >::iterator it = imagePoints.begin(); it != imagePoints.end(); ++it) {

    //     int cnt2 = 0;
    //     for (vector<Point2f>::iterator it2 = (*it).begin(); it2 != (*it).end(); ++it2) {
    //         printf("  %d-%d: %f, %f\n", cnt1, cnt2, it2->x, it2->y);
    //         cnt2 ++;
    //     }
    //     cnt1 ++;
    // }
    // printf("---------------------\n");
    // // vector<vector<Point3f> > objectPoints
    // cnt1 = 0;
    // for (vector<vector<Point3f> >::iterator it = objectPoints.begin(); it != objectPoints.end(); ++it) {

    //     int cnt2 = 0;
    //     for (vector<Point3f>::iterator it2 = (*it).begin(); it2 != (*it).end(); ++it2) {
    //         printf("  %d-%d: %f, %f\n", cnt1, cnt2, it2->x, it2->y);
    //         cnt2 ++;
    //     }
    //     cnt1 ++;
    // }



    double rms = calibrateCamera(objectPoints, imagePoints, imageSize, cameraMatrix,
                    distCoeffs, rvecs, tvecs, flags|CALIB_FIX_K4|CALIB_FIX_K5);
                    ///*|CALIB_FIX_K3*/|CALIB_FIX_K4|CALIB_FIX_K5);
    printf("RMS error reported by calibrateCamera: %g\n", rms);

    bool ok = checkRange(cameraMatrix) && checkRange(distCoeffs);

    totalAvgErr = computeReprojectionErrors(objectPoints, imagePoints,
                rvecs, tvecs, cameraMatrix, distCoeffs, reprojErrs);


    return ok;
}


static void saveCameraParams( const string& filename,
                       Size imageSize, Size boardSize,
                       float squareSize, float aspectRatio, int flags,
                       const Mat& cameraMatrix, const Mat& distCoeffs,
                       const vector<Mat>& rvecs, const vector<Mat>& tvecs,
                       const vector<float>& reprojErrs,
                       const vector<vector<Point2f> >& imagePoints,
                       double totalAvgErr )
{
    FileStorage fs( filename, FileStorage::WRITE );

    time_t tt;
    time( &tt );
    struct tm *t2 = localtime( &tt );
    char buf[1024];
    strftime( buf, sizeof(buf)-1, "%c", t2 );

    fs << "calibration_time" << buf;

    if( !rvecs.empty() || !reprojErrs.empty() )
        fs << "nframes" << (int)std::max(rvecs.size(), reprojErrs.size());
    fs << "image_width" << imageSize.width;
    fs << "image_height" << imageSize.height;
    fs << "board_width" << boardSize.width;
    fs << "board_height" << boardSize.height;
    fs << "square_size" << squareSize;

    if( flags & CALIB_FIX_ASPECT_RATIO )
        fs << "aspectRatio" << aspectRatio;

    if( flags != 0 )
    {
        sprintf( buf, "flags: %s%s%s%s",
            flags & CALIB_USE_INTRINSIC_GUESS ? "+use_intrinsic_guess" : "",
            flags & CALIB_FIX_ASPECT_RATIO ? "+fix_aspectRatio" : "",
            flags & CALIB_FIX_PRINCIPAL_POINT ? "+fix_principal_point" : "",
            flags & CALIB_ZERO_TANGENT_DIST ? "+zero_tangent_dist" : "" );
        //cvWriteComment( *fs, buf, 0 );
    }

    fs << "flags" << flags;

    fs << "camera_matrix" << cameraMatrix;
    fs << "distortion_coefficients" << distCoeffs;

    fs << "avg_reprojection_error" << totalAvgErr;
    if( !reprojErrs.empty() )
        fs << "per_view_reprojection_errors" << Mat(reprojErrs);

    if( !rvecs.empty() && !tvecs.empty() )
    {
        CV_Assert(rvecs[0].type() == tvecs[0].type());
        Mat bigmat((int)rvecs.size(), 6, rvecs[0].type());
        for( int i = 0; i < (int)rvecs.size(); i++ )
        {
            Mat r = bigmat(Range(i, i+1), Range(0,3));
            Mat t = bigmat(Range(i, i+1), Range(3,6));

            CV_Assert(rvecs[i].rows == 3 && rvecs[i].cols == 1);
            CV_Assert(tvecs[i].rows == 3 && tvecs[i].cols == 1);
            //*.t() is MatExpr (not Mat) so we can use assignment operator
            r = rvecs[i].t();
            t = tvecs[i].t();
        }
        //cvWriteComment( *fs, "a set of 6-tuples (rotation vector + translation vector) for each view", 0 );
        fs << "extrinsic_parameters" << bigmat;
    }

    if( !imagePoints.empty() )
    {
        Mat imagePtMat((int)imagePoints.size(), (int)imagePoints[0].size(), CV_32FC2);
        for( int i = 0; i < (int)imagePoints.size(); i++ )
        {
            Mat r = imagePtMat.row(i).reshape(2, imagePtMat.cols);
            Mat imgpti(imagePoints[i]);
            imgpti.copyTo(r);
        }
        fs << "image_points" << imagePtMat;
    }
}

static bool readStringList( const string& filename, vector<string>& l )
{
    l.resize(0);
    FileStorage fs(filename, FileStorage::READ);
    if( !fs.isOpened() )
        return false;
    FileNode n = fs.getFirstTopLevelNode();
    if( n.type() != FileNode::SEQ )
        return false;
    FileNodeIterator it = n.begin(), it_end = n.end();
    for( ; it != it_end; ++it )
        l.push_back((string)*it);
    return true;
}


static bool runAndSave(const string& outputFilename,
                const vector<vector<Point2f> >& imagePoints,
                Size imageSize, Size boardSize, Pattern patternType, float squareSize,
                float aspectRatio, int flags, Mat& cameraMatrix,
                Mat& distCoeffs, 
                vector<Mat>& rvecs, vector<Mat>& tvecs,
                bool writeExtrinsics, bool writePoints )
{
    // vector<Mat> rvecs, tvecs;
    vector<float> reprojErrs;
    double totalAvgErr = 0;

    bool ok = runCalibration(imagePoints, imageSize, boardSize, patternType, squareSize,
                   aspectRatio, flags, cameraMatrix, distCoeffs,
                   rvecs, tvecs, reprojErrs, totalAvgErr);
    printf("%s. avg reprojection error = %.2f\n",
           ok ? "Calibration succeeded" : "Calibration failed",
           totalAvgErr);

    if( ok )
        saveCameraParams( outputFilename, imageSize,
                         boardSize, squareSize, aspectRatio,
                         flags, cameraMatrix, distCoeffs,
                         writeExtrinsics ? rvecs : vector<Mat>(),
                         writeExtrinsics ? tvecs : vector<Mat>(),
                         writeExtrinsics ? reprojErrs : vector<float>(),
                         writePoints ? imagePoints : vector<vector<Point2f> >(),
                         totalAvgErr );
    return ok;
}


// ycc
static void projectPixel2World(Mat trans_mat, Point2f image_point, Point2f& world_point)
{

    double x1 = trans_mat.at<double>(0,0) - trans_mat.at<double>(2,0) * image_point.x;
    double y1 = trans_mat.at<double>(0,1) - trans_mat.at<double>(2,1) * image_point.x;
    double c1 = trans_mat.at<double>(2,3) * image_point.x - trans_mat.at<double>(0,3);
    double x2 = trans_mat.at<double>(1,0) - trans_mat.at<double>(2,0) * image_point.y;
    double y2 = trans_mat.at<double>(1,1) - trans_mat.at<double>(2,1) * image_point.y;
    double c2 = trans_mat.at<double>(2,3) * image_point.y - trans_mat.at<double>(1,3);
    world_point.y = (c1*x2 - c2*x1) / (x2*y1 - x1*y2);
    world_point.x = (c1 - y1*world_point.y) / x1;
}


// ycc
void CalibrationIlluminane(Mat &image, double ratio)
{
    image = image * ratio;
}


int main( int argc, char** argv )
{
    Size boardSize, imageSize;
    float squareSize = 1.f, aspectRatio = 1.f;
    Mat cameraMatrix, distCoeffs;
    const char* outputFilename = "out_camera_data.yml";
    const char* inputFilename = 0;

    int i, nframes = 2;
    bool writeExtrinsics = false, writePoints = false;
    bool undistortImage = false;
    int flags = 0;
    VideoCapture capture;
    bool flipVertical = false;
    bool showUndistorted = false;
    bool videofile = false;
    int delay = 1000;
    clock_t prevTimestamp = 0;
    int mode = DETECTION;
    int cameraId = 0;
    vector<vector<Point2f> > imagePoints;
    vector<string> imageList;
    Pattern pattern = CHESSBOARD;

    // ycc
    vector<Mat> rvecs, tvecs;
    vector<vector<Point3f> > objectPoints(1);
    Mat rmat, rtmat, trans_mat;


    if( argc < 2 )
    {
        help();
        return 0;
    }

    for( i = 1; i < argc; i++ )
    {
        const char* s = argv[i];
        if( strcmp( s, "-w" ) == 0 )
        {
            if( sscanf( argv[++i], "%u", &boardSize.width ) != 1 || boardSize.width <= 0 )
                return fprintf( stderr, "Invalid board width\n" ), -1;
        }
        else if( strcmp( s, "-h" ) == 0 )
        {
            if( sscanf( argv[++i], "%u", &boardSize.height ) != 1 || boardSize.height <= 0 )
                return fprintf( stderr, "Invalid board height\n" ), -1;
        }
        else if( strcmp( s, "-pt" ) == 0 )
        {
            i++;
            if( !strcmp( argv[i], "circles" ) )
                pattern = CIRCLES_GRID;
            else if( !strcmp( argv[i], "acircles" ) )
                pattern = ASYMMETRIC_CIRCLES_GRID;
            else if( !strcmp( argv[i], "chessboard" ) )
                pattern = CHESSBOARD;
            else
                return fprintf( stderr, "Invalid pattern type: must be chessboard or circles\n" ), -1;
        }
        else if( strcmp( s, "-s" ) == 0 )
        {
            if( sscanf( argv[++i], "%f", &squareSize ) != 1 || squareSize <= 0 )
                return fprintf( stderr, "Invalid board square width\n" ), -1;
        }
        else if( strcmp( s, "-n" ) == 0 )
        {
            if( sscanf( argv[++i], "%u", &nframes ) != 1 || nframes <= 3 )
                return printf("Invalid number of images\n" ), -1;
        }
        else if( strcmp( s, "-a" ) == 0 )
        {
            if( sscanf( argv[++i], "%f", &aspectRatio ) != 1 || aspectRatio <= 0 )
                return printf("Invalid aspect ratio\n" ), -1;
            flags |= CALIB_FIX_ASPECT_RATIO;
        }
        else if( strcmp( s, "-d" ) == 0 )
        {
            if( sscanf( argv[++i], "%u", &delay ) != 1 || delay <= 0 )
                return printf("Invalid delay\n" ), -1;
        }
        else if( strcmp( s, "-op" ) == 0 )
        {
            writePoints = true;
        }
        else if( strcmp( s, "-oe" ) == 0 )
        {
            writeExtrinsics = true;
        }
        else if( strcmp( s, "-zt" ) == 0 )
        {
            flags |= CALIB_ZERO_TANGENT_DIST;
        }
        else if( strcmp( s, "-p" ) == 0 )
        {
            flags |= CALIB_FIX_PRINCIPAL_POINT;
        }
        else if( strcmp( s, "-v" ) == 0 )
        {
            flipVertical = true;
        }
        else if( strcmp( s, "-V" ) == 0 )
        {
            videofile = true;
        }
        else if( strcmp( s, "-o" ) == 0 )
        {
            outputFilename = argv[++i];
        }
        else if( strcmp( s, "-su" ) == 0 )
        {
            showUndistorted = true;
        }
        else if( s[0] != '-' )
        {
            if( isdigit(s[0]) )
                sscanf(s, "%d", &cameraId);
            else
                inputFilename = s;
        }
        else
            return fprintf( stderr, "Unknown option %s", s ), -1;
    }

    if( inputFilename )
    {
        if( !videofile && readStringList(inputFilename, imageList) )
            mode = CAPTURING;
        else
            capture.open(inputFilename);
    }
    else
        capture.open(cameraId);

    if( !capture.isOpened() && imageList.empty() )
        return fprintf( stderr, "Could not initialize video (%d) capture\n",cameraId ), -2;

    if( !imageList.empty() )
        nframes = (int)imageList.size();

    if( capture.isOpened() )
        printf( "%s", liveCaptureHelp );

    namedWindow( "Image View");

    for(i = 0;;i++)
    {
        Mat view, viewGray;
        bool blink = false;

        if( capture.isOpened() )
        {
            Mat view0;
            capture >> view0;
            view0.copyTo(view);
        }
        else if( i < (int)imageList.size() )
            view = imread(imageList[i], 1);

        if(view.empty())
        {
            if( imagePoints.size() > 0 )
                runAndSave(outputFilename, imagePoints, imageSize,
                           boardSize, pattern, squareSize, aspectRatio,
                           flags, cameraMatrix, distCoeffs,
                           rvecs, tvecs,
                           writeExtrinsics, writePoints);
            break;
        }

        imageSize = view.size();

        if( flipVertical )
            flip( view, view, 0 );

        vector<Point2f> pointbuf;
        cvtColor(view, viewGray, COLOR_BGR2GRAY);
        

        // ycc
        // view = viewGray;
        // double alpha = 1;
        // double beta = 100;
        // view = alpha * view + beta;


        bool found;
        switch( pattern )
        {
            case CHESSBOARD:
                found = findChessboardCorners( view, boardSize, pointbuf,
                    CALIB_CB_ADAPTIVE_THRESH | CALIB_CB_FAST_CHECK | CALIB_CB_NORMALIZE_IMAGE);
                break;
            case CIRCLES_GRID:
                found = findCirclesGrid( view, boardSize, pointbuf );
                break;
            case ASYMMETRIC_CIRCLES_GRID:
                found = findCirclesGrid( view, boardSize, pointbuf, CALIB_CB_ASYMMETRIC_GRID );
                break;
            default:
                return fprintf( stderr, "Unknown pattern type\n" ), -1;
        }

       // improve the found corners' coordinate accuracy
        if( pattern == CHESSBOARD && found) cornerSubPix( viewGray, pointbuf, Size(11,11),
            Size(-1,-1), TermCriteria( TermCriteria::EPS+TermCriteria::COUNT, 30, 0.1 ));

        if( mode == CAPTURING && found &&
           (!capture.isOpened() || clock() - prevTimestamp > delay*1e-3*CLOCKS_PER_SEC) )
        {
            imagePoints.push_back(pointbuf);
            prevTimestamp = clock();
            blink = capture.isOpened();
        }

        if(found) {
            drawChessboardCorners( view, boardSize, Mat(pointbuf), found );

            // ycc
            // Label coordinate
            if(mode == CAPTURING || mode == DETECTION) {
                int cnt = 0;
                for(int i = 0; i < pointbuf.size(); i ++) {
                    char tmpx[10];
                    char tmpy[10];
                    sprintf(tmpx, "%.2f", pointbuf[i].x);
                    sprintf(tmpy, "%.2f", pointbuf[i].y);
                    Point textOrigin(pointbuf[i].x, pointbuf[i].y);
                    String msg = "(" + String(tmpx) + "," + String(tmpy) + ")";
                    putText( view, msg, textOrigin, 1, 1,
                        mode != CALIBRATED ? Scalar(0,0,255) : Scalar(0,255,0));
                }
            }

        }

        string msg = mode == CAPTURING ? "100/100" :
            mode == CALIBRATED ? "Calibrated" : "Press 'g' to start";
        int baseLine = 0;
        Size textSize = getTextSize(msg, 1, 1, 1, &baseLine);
        Point textOrigin(view.cols - 2*textSize.width - 10, view.rows - 2*baseLine - 10);
        
        if( mode == CAPTURING )
        {
            if(undistortImage)
                msg = format( "%d/%d Undist", (int)imagePoints.size(), nframes );
            else
                msg = format( "%d/%d", (int)imagePoints.size(), nframes );
        }

        putText( view, msg, textOrigin, 1, 1,
                 mode != CALIBRATED ? Scalar(0,0,255) : Scalar(0,255,0));

        if( blink )
            bitwise_not(view, view);

        if( mode == CALIBRATED && undistortImage )
        {
            Mat temp = view.clone();
            undistort(temp, view, cameraMatrix, distCoeffs);
        }

        
        // ycc
        if(found && mode == CALIBRATED) {
            int cnt = 0;

            // printf("================\n");
            // printf("rvecs\n");
            // cnt = 0;
            // for(int i = 0; i < rvecs.size(); i ++) {
            //     printf("--------------\n");
            //     printf("  %d: rows=%d, cols=%d\n", cnt, rvecs[i].rows, rvecs[i].cols);
            //     cout << "  rvecs = " << endl << rvecs[i] << endl;
            //     cnt ++;
            // }
            // printf("================\n");
            // printf("tvecs\n");
            // cnt = 0;
            // for(int i = 0; i < tvecs.size(); i ++) {
            //     printf("--------------\n");
            //     printf("  %d: rows=%d, cols=%d\n", cnt, tvecs[i].rows, tvecs[i].cols);
            //     cout << "  tvecs = " << endl << tvecs[i] << endl;
            //     cnt ++;
            // }

            // printf("================\n");
            // printf("coners\n");
            

            // use projectPoints to map from object to image coordinate
            // vector<Point2f> imagePoints2;
            // projectPoints(Mat(objectPoints[0]), rvecs[0], tvecs[0],
            //       cameraMatrix, distCoeffs, imagePoints2);


            // use m = A[R|t]M to map from object to image coordinate
            // vector<Point2f> imagePoints3;
            vector<Point2f> objectPoints3;
            for(int i = 0; i < objectPoints[0].size(); i ++) {
                
                // map
                // Mat ret1;
                // Mat tmp_obj(4, 1, cameraMatrix.type());
                // tmp_obj.at<double>(0,0) = objectPoints[0][i].x;
                // tmp_obj.at<double>(1,0) = objectPoints[0][i].y;
                // tmp_obj.at<double>(2,0) = 0;
                // tmp_obj.at<double>(3,0) = 1;

                // ret1 = trans_mat * tmp_obj;
                // Point2f tmp2(ret1.at<double>(0,0)/ret1.at<double>(2,0), ret1.at<double>(1,0)/ret1.at<double>(2,0));
                // imagePoints3.push_back(tmp2);


                // map back
                Mat ret2;
                Mat tmp_img(3, 1, cameraMatrix.type());
                // tmp_img.at<double>(0,0) = imagePoints[0][i].x;
                // tmp_img.at<double>(1,0) = imagePoints[0][i].y;
                // tmp_img.at<double>(2,0) = 1;

                // double x1 = trans_mat.at<double>(0,0) - trans_mat.at<double>(2,0) * tmp_img.at<double>(0,0);
                // double y1 = trans_mat.at<double>(0,1) - trans_mat.at<double>(2,1) * tmp_img.at<double>(0,0);
                // double c1 = trans_mat.at<double>(2,3) * tmp_img.at<double>(0,0) - trans_mat.at<double>(0,3);
                // double x2 = trans_mat.at<double>(1,0) - trans_mat.at<double>(2,0) * tmp_img.at<double>(1,0);
                // double y2 = trans_mat.at<double>(1,1) - trans_mat.at<double>(2,1) * tmp_img.at<double>(1,0);
                // double c2 = trans_mat.at<double>(2,3) * tmp_img.at<double>(1,0) - trans_mat.at<double>(1,3);
                // double Y = (c1*x2 - c2*x1) / (x2*y1 - x1*y2);
                // double X = (c1 - y1*Y) / x1;
                // double X2 = (c2 - y2*Y) / x2;
                Point2f tmp3;
                projectPixel2World(trans_mat, imagePoints[0][i], tmp3);
                // printf("X=%f, Y=%f\n", tmp3.x, tmp3.y);
                
                objectPoints3.push_back(tmp3);
            }

            
            cnt = 0;
            for(int i = 0; i < pointbuf.size(); i ++) {
                char tmpx[10];
                char tmpy[10];
                sprintf(tmpx, "%.2f", pointbuf[i].x);
                sprintf(tmpy, "%.2f", pointbuf[i].y);
                String msg1 = "(" + String(tmpx) + "," + String(tmpy) + ")";
                
                Point textOrigin1(pointbuf[i].x, pointbuf[i].y);
                putText( view, msg1, textOrigin1, 1, 1, Scalar(255,0,0));


                char tmpX[10];
                char tmpY[10];
                sprintf(tmpX, "%.2f", objectPoints[0][i].x);
                sprintf(tmpY, "%.2f", objectPoints[0][i].y);
                String msg2 = "(" + String(tmpX) + "," + String(tmpY) + ")";

                int baseLine;
                Size textSize1 = getTextSize(msg1, 1, 1, 1, &baseLine);
                Point textOrigin2(pointbuf[i].x, pointbuf[i].y - 1.5*textSize1.height);
                putText( view, msg2, textOrigin2, 1, 1, Scalar(255,0,0));

                // printf("  %d: (%.2f,%.2f) - (%.2f,%.2f) - (%.2f,%.2f) - (%.2f,%.2f) - (%.2f,%.2f)\n", i,
                //     objectPoints[0][i].x, objectPoints[0][i].y,
                //     pointbuf[i].x, pointbuf[i].y,
                //     imagePoints2[i].x, imagePoints2[i].y,
                //     imagePoints3[i].x, imagePoints3[i].y,
                //     objectPoints3[i].x, objectPoints3[i].y);
                printf("  %d: (%.2f,%.2f) - (%.2f,%.2f) - (%.2f,%.2f)\n", i,
                    objectPoints[0][i].x, objectPoints[0][i].y,
                    pointbuf[i].x, pointbuf[i].y,
                    objectPoints3[i].x, objectPoints3[i].y);
            }

        }

        if(mode == CALIBRATED) {

            // tracking
            Mat viewHSV;
            cvtColor(view, viewHSV, COLOR_BGR2HSV);
            // define range of blue color in    
            Scalar lower_blue(110,50,50);
            Scalar upper_blue(130,255,255);
            Mat mask;
            inRange(viewHSV, lower_blue, upper_blue, mask);
            // view = mask;
            
            Moments mu = moments(mask);
            Point2f ctr(mu.m10/mu.m00 , mu.m01/mu.m00 );
            char tmp[100];
            sprintf(tmp, "(%.2f,%.2f)\n", ctr.x, ctr.y);
            String ctr_msg = String(tmp);
            putText(view, ctr_msg, ctr, 1, 1, Scalar(0,0,255));

            Point2f ctr_world;
            projectPixel2World(trans_mat, ctr, ctr_world);
            sprintf(tmp, "(%.2f,%.2f)\n", ctr_world.x, ctr_world.y);
            String ctr_world_msg = String(tmp);
            int baseLine;
            Size textSize1 = getTextSize(ctr_msg, 1, 1, 1, &baseLine);
            Point textOrigin2(ctr.x, ctr.y - 1.5*textSize1.height);
            putText(view, ctr_world_msg, textOrigin2, 1, 1, Scalar(0,0,255));
            
        }


        imshow("Image View", view);
        // imshow("Image View", mask);
        // imshow("Image View", res);


        // ycc
        if(found && mode == CALIBRATED) {
            printf(">> enter to continue\n");
            waitKey(0);
        }


        int key = 0xff & waitKey(capture.isOpened() ? 50 : 500);

        if( (key & 255) == 27 )
            break;

        if( key == 'u' && mode == CALIBRATED )
            undistortImage = !undistortImage;

        if( capture.isOpened() && key == 'g' )
        {
            mode = CAPTURING;
            imagePoints.clear();
        }

        if( mode == CAPTURING && imagePoints.size() >= (unsigned)nframes )
        {
            if( runAndSave(outputFilename, imagePoints, imageSize,
                       boardSize, pattern, squareSize, aspectRatio,
                       flags, cameraMatrix, distCoeffs,
                       rvecs, tvecs,
                       writeExtrinsics, writePoints)) {

                mode = CALIBRATED;

                // ycc
                calcChessboardCorners(boardSize, squareSize, objectPoints[0]);
                objectPoints.resize(imagePoints.size(),objectPoints[0]);
                // std::cout << "  cameraMatrix = " << endl << cameraMatrix << endl;

                Rodrigues(Mat(rvecs[0]), rmat);
                // std::cout << "r mat = " << endl << rmat << endl;

                hconcat(rmat, tvecs[0], rtmat);
                // std::cout << "rt mat = " << endl << rtmat << endl;

                cameraMatrix.copyTo(trans_mat);
                // cout << "trans mat = " << endl << trans_mat << endl;
                trans_mat = trans_mat * rtmat;
                // cout << "trans mat = " << endl << trans_mat << endl;
            }
            else
                mode = DETECTION;
            if( !capture.isOpened() )
                break;
        }
    }


    if( !capture.isOpened() && showUndistorted )
    {
        Mat view, rview, map1, map2;
        initUndistortRectifyMap(cameraMatrix, distCoeffs, Mat(),
                                getOptimalNewCameraMatrix(cameraMatrix, distCoeffs, imageSize, 1, imageSize, 0),
                                imageSize, CV_16SC2, map1, map2);

        for( i = 0; i < (int)imageList.size(); i++ )
        {
            view = imread(imageList[i], 1);
            if(view.empty())
                continue;
            //undistort( view, rview, cameraMatrix, distCoeffs, cameraMatrix );
            remap(view, rview, map1, map2, INTER_LINEAR);
            imshow("Image View", rview);
            int c = 0xff & waitKey();
            if( (c & 255) == 27 || c == 'q' || c == 'Q' )
                break;
        }
    }

    return 0;
}
