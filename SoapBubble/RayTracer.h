//
//  RayTracer.h
//  SoapBubble
//
//  Created by Kevin Utecht on 6/21/13.
//
//  Modeling Interference Color
//
//  Created for Computer Graphics course CS 5117 in 1993.  Updated in 2013 to turn into an iPad app.
//
//  The ray tracing algorithm is derived from code by George Kyriazis ((c) 1988 by George Kyriazis)
//  ( see http://cd.textfiles.com/graphics16000/GENERAL/RAYTRACE/RAY_20/RAY2_0.TXT )
//
//  The fresnel_formula was described in article Ray Tracing Interference Color by Maria Lurdes Dias in
//  1991 magazine IEEE Computer Graphics & Applications.  Maria had a follow-up article
//  Ray Tracing Interference Color: Visualizing Newton's Ringswith some additional clarification
//  in 1994.
//

#import <UIKit/UIKit.h>


UIImage* raytrace(int xRes, int yRes, double bubbleThicknessCoefficient);