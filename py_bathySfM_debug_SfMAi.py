# -*- coding: utf-8 -*-
"""
Created on Wed Nov  7 15:47:38 2018

@author: jdietric
"""

#pyqt import
#from PyQt5 import QtCore, QtGui, QtWidgets

# other imports
import os
import sys
import numpy as np
import pandas as pd
import sympy.geometry as spg
import matplotlib.path as mplPath
from datetime import datetime

# MAIN PROGRAM HELPER FUNCTIONS (run on OK button)

def footprints(cam, sensor, base_elev):
    """
    This function calculates the instantaneous field of view (IFOV) for 
    the camera(s) that are passed.\n
    Vars:\n
    \t cam = pandas dataframe (n x ~6, fields: x,y,z,yaw,pitch,roll)\n
    \t sensor = pandas dataframe (1 x 3, fields: focal, sensor_x, sensor_y):
    \t focal length (mm), sensor x dim (mm), sensor y dim (mm)\n
    \t base_elev = average elevation of your site (meters, or in the same
    \t measure as your coordinates)\n
    Creates approx. coordinates for sensor
    corners (north-oriented and zero pitch) at the camera's x,y,z. Rotates
    the sensor coords in 3D space to the camera's pitch and yaw angles (roll
    angles are ignored for now) and projects corner rays through the camera 
    x,y,z to a approx ground plane. The intersection of the rays with the
    ground are the corners of the photo footprint.\n
    *** Photos that have picth angles that cause the horizon to be visable will
    cause the UL and UR path coordniates to wrong. These cameras are 
    disreguarded and the footprint will be set to NaN in the output.***\n 
    RETURNS: footprints = Pandas dataframe (n x 1) of Matplotlib Path objects()
    """
    
    # Setup DF to house camera footprint polygons
    footprints = pd.DataFrame(np.zeros((cam.shape[0],1)), columns=['fov'])
    
    # debug - blank 3d array for inter_points
#        itp_f = '//thor.ad.uni.edu/users/jdietric/Documents/Python Scripts/py_sfm_depth/WhiteR_2016/itp.npy'
#        itp = np.zeros((cam.shape[0],4,2))
    
    # convert sensor dimensions to meters, divide x/y for corner coord calc
    f = sensor.focal[0] * 0.001
    sx = sensor.sensor_x[0] / 2 * 0.001
    sy = sensor.sensor_y[0] / 2 * 0.001

    # calculate the critical pitch (in degrees) where the horizon will be 
    #   visible with the horizon viable, the ray projections go backward 
    #   and produce erroneous IFOV polygons (90 - 0.5*vert_fov)
    crit_pitch = 90 - np.rad2deg(np.arctan(sy / f))
    
    # User Feedback
    print("Proccesing Camera IFOVs (%i total)..." %(cam.shape[0]))
    sys.stdout.flush()
     
    # for each camera...
    for idx, row in cam.iterrows():
        
        # check is the camera pitch is over the critical value
        if row.pitch < crit_pitch:
            
            # sensor corners (UR,LR,LL,UL), north-oriented and zero pitch
            corners = np.array([[row.x+sx,row.y-f,row.z+sy],
                               [row.x+sx,row.y-f,row.z-sy],
                               [row.x-sx,row.y-f,row.z-sy],
                               [row.x-sx,row.y-f,row.z+sy]])
            
            # offset corner points by cam x,y,z for rotation
            cam_pt = np.atleast_2d(np.array([row.x, row.y, row.z]))
            corner_p = corners - cam_pt
    
            # get pitch and yaw from the camera, convert to radians
            pitch = np.deg2rad(90.0-row.pitch)
            roll = np.deg2rad(row.roll)
            yaw = np.deg2rad(row.yaw)
            
            # setup picth rotation matrix (r_x) and yaw rotation matrix (r_z)
            r_x = np.matrix([[1.0,0.0,0.0],
                             [0.0,np.cos(pitch),-1*np.sin(pitch)],
                             [0.0,np.sin(pitch),np.cos(pitch)]])
                             
            r_y = np.matrix([[np.cos(roll),0.0,np.sin(roll)],
                             [0.0,1.0,0.0],
                             [-1*np.sin(roll),0.0,np.cos(roll)]])
            
            r_z =  np.matrix([[np.cos(yaw),-1*np.sin(yaw),0],
                              [np.sin(yaw),np.cos(yaw),0],
                              [0,0,1]])
            
            # rotate corner_p by r_x, then r_z, add back cam x,y,z offsets
            # produces corner coords rotated for pitch and yaw
            p_pr = np.matmul(np.matmul(corner_p, r_x),r_y)            
            p_out = np.matmul(p_pr, r_z) + cam_pt
            
            # GEOMETRY
            # Set Sympy 3D point for the camera and a 3D plane for intersection
            cam_sp = spg.Point3D(row.x, row.y, row.z)
            plane = spg.Plane(spg.Point3D(row.x, row.y, base_elev),
                                      normal_vector=(0,0,1))
            
            # blank array for footprint intersection coords
            inter_points = np.zeros((corners.shape[0],2))
            
            # for each sensor corner point
            idx_b = 0
            for pt in np.asarray(p_out):
                
                # create a Sympy 3D point and create a Sympy 3D ray from 
                #   corner point through camera point
                pt_sp = spg.Point3D(pt[0],pt[1],pt[2])
                ray = spg.Ray3D(pt_sp,cam_sp)
                
                # calculate the intersection of the ray with the plane                
                inter_pt = plane.intersection(ray)
                
                # Extract out the X,Y coords fot eh intersection point
                #   ground intersect points will be in this order (LL,UL,UR,LR)
                inter_points[idx_b,0] = inter_pt[0].x.evalf()
                inter_points[idx_b,1] = inter_pt[0].y.evalf()
                
                idx_b += 1
        
        # if crit_pitch is exceeded set inter_points to NaN
        else:
            inter_points = np.full((4,2),np.nan)
        
        # append inter_points to footprints as a matplotlib path object
        footprints.fov[idx] = mplPath.Path(inter_points)
        
        #debug - save inter_points
#            itp[idx,:,:] = inter_points
        
        # User feedback and progress bar
        if (idx+1) % 10 == 0:
            print("%i cameras processed..." %(idx+1))
#            gui.top_progBar.setValue(idx)
#            gui.topProg_Lbl.setText("Calculating Camera Footprints - " + str(idx+1))
            #sys.stdout.flush()
            
    #debug - save inter_points
    #np.save(itp_f,itp)
    
    return footprints
# END - footprints
    
def visibility(cam, sensor, footprints, targets):
    """    
    This function tests is the target points (x,y only) are "visable" (i.e.
    within the photo footprints) and calculates the "r" angle for the refraction 
    correction\n
    Vars:\n
    \t cam = Pandas dataframe (n x ~6, fields: x,y,z,yaw,pitch,roll)\n
    \t footprints = Pandas dataframe (n x 1) of Matplotlib Path objects\n
    \t targets = Pandas dataframe (n x ~3, fields: x,y,sfm_z...)\n
    
    RETURNS: r_filt = numpy array (n_points x n_cams) of filtered "r" angles.\n
    Points that are not visable to a camera will have a NaN "r" angle. 
    """
    
    # Setup boolean array for visability
    vis = np.zeros((targets.shape[0],cam.shape[0])) 
    
    # for each path objec in footprints, check is the points in targets are
    #   within the path polygon. path.contains_points returns boolean.
    #   the results are accumulated in the vis array.
    for idx in range(footprints.shape[0]):
        path = footprints.fov[idx]
        vis[:,idx] = path.contains_points(np.array([targets.x.values, targets.y.values]).T)
    
    # calculate the coord. deltas between the cameras and the target
    dx = np.atleast_2d(cam.x.values) - np.atleast_2d(targets.x.values).T
    dy = np.atleast_2d(cam.y.values) - np.atleast_2d(targets.y.values).T
    dz = np.atleast_2d(cam.z.values) - np.atleast_2d(targets.z.values).T
    
    # calc xy distance (d)
    d = np.sqrt((dx)**2+(dy)**2)
    
    # calc inclination angle (r) from targets to cams
    r = np.rad2deg(np.arctan(d/dz))
    
    # radial distance from the image center
    # calc mm/pix in the image size, convert focal length to pixels
    mm_pix = ((sensor.sensor_x/sensor.pix_x)+(sensor.sensor_y/sensor.pix_y))/2
    f_pix = sensor.focal / mm_pix[0]
    rad_dist = f_pix[0] * np.tan(np.radians(r))
    rad_dist = rad_dist * vis
    rad_dist[rad_dist == 0] = np.nan
    
    max_radial_dist = sensor.pix_x[0]
    rad_dist_pcent = rad_dist / max_radial_dist
    
    # slope dist calc.
    # cosine
    SD = dz/np.cos(np.radians(r))
    
    #filter r, d, and sd values by vis matrix
    r_filt = r * vis
    r_filt[r_filt == 0] = np.nan
    
    d_filt = d * vis
    d_filt[d_filt == 0] = np.nan
    
    SD_filt = SD * vis
    SD_filt[SD_filt == 0] = np.nan
    
    if 'quality' in cam.columns:
        qual = cam.quality.values
        qual_mat = np.repeat(np.atleast_2d(qual),targets.shape[0],axis=0)
        qual_vis = qual_mat * vis
        qual_vis[qual_vis==0] = np.nan
    else:      
        qual_vis = np.zeros((targets.shape[0],cam.shape[0]))
        qual_vis[:] = np.nan
        
    return r_filt, d_filt, SD_filt, rad_dist, qual_vis 

#def correction(r, target, extras):
#    """Performs the per camera refraction correction on a target point.
#    Refer to the documentation for the specifics.
#    extras = [smAng, stats, camStats, weight]"""
#    
#    # convert r array to radians for trig calculations
#    ang_r = np.radians(r)
#
#    # calculate the refraction angle i 
#    ang_i = np.arcsin(1.0/1.337 * np.sin(ang_r))
#    
#    # calculate the apparent depth from the water surface elev. and the 
#    #   target SfM elevation
#    target['h_a'] = target.w_surf - target.sfm_z
#
#    # calculate the distance from the point to the air/water interface
#    x_dist = np.array([target.h_a.values]).T * np.tan(ang_r)
#    
#    # calculate the corrected (actual) depth
#    h =  x_dist / np.tan(ang_i)
#   
#    # subtract the corrected depth from the water surface elevation to get the
#    #   corrected elevation
#    cor_elev = np.array([target.w_surf]).T - h
#    
#    # append the mean values for the actual depth and corrected elevation to
#    #   the target data frame   
#    target['h_avg'] = np.nanmean(h, axis = 1)
#    target['corElev_avg'] = np.nanmean(cor_elev, axis = 1)
#    
#    if extras[0]:
#    # calc small angle approximation
#        target['smAng_h'] = target['h_a'] * 1.34
#        target['smAng_elev'] = target.w_surf - target['smAng_h']       
#       
    #if extras[1]:
    # some extra statistics for playing around (option from GUI)
#        target['h_std'] = np.nanstd(h, axis = 1)
#        target['h_med'] = np.nanmedian(h, axis = 1)
#        target['h_min'] = np.nanmin(h, axis = 1)
#        target['h_max'] = np.nanmax(h, axis = 1)
#        target['h_per15'] = np.nanpercentile(h, 15,axis = 1)
#        target['h_per25'] = np.nanpercentile(h, 25,axis = 1)
#        target['h_per55'] = np.nanpercentile(h, 55,axis = 1)
#        target['h_per60'] = np.nanpercentile(h, 60,axis = 1)
#        target['h_per75'] = np.nanpercentile(h, 75,axis = 1)
#        target['h_per80'] = np.nanpercentile(h, 80,axis = 1)
#        target['h_per90'] = np.nanpercentile(h, 90,axis = 1)
#        target['h_per95'] = np.nanpercentile(h, 95,axis = 1)
#        target['h_iqr'] = target['h_per75'] - target['h_per25']
#        target['h_lif'] = target['h_per25'] - (1.5 * target['h_iqr'])
#        target['h_uif'] = target['h_per75'] + (1.5 * target['h_iqr'])
#        target['h_lof'] = target['h_per25'] - (3 * target['h_iqr'])
#        target['h_uof'] = target['h_per75'] + (3 * target['h_iqr'])
#        l_whisk = np.repeat(np.atleast_2d(target['h_lif']).T,188,axis=1)
#        u_whisk = np.repeat(np.atleast_2d(target['h_uif']).T,188,axis=1)
#        target['u_whisk'] = np.nanmax(np.ma.masked_less_equal(h,u_whisk).data,axis=1)
#        target['l_whisk'] = np.nanmin(np.ma.masked_greater_equal(h,l_whisk).data,axis=1)
#        mild_out = np.zeros_like(h)
#        ext_out = np.zeros_like(h)
#        mild_out[h < np.atleast_2d(target['h_lif']).T] = 1
#        mild_out[h > np.atleast_2d(target['h_uif']).T] = 1   
#        ext_out[h < np.atleast_2d(target['h_lof']).T] = 1
#        ext_out[h > np.atleast_2d(target['h_uof']).T] = 1   
#        target['h_mildout'] = np.nansum(mild_out, axis = 1)
#        target['h_extout'] = np.nansum(ext_out, axis = 1)
        
#    if extras[3]:
#    # if the weighted checkbox is selected
#        max_angle = np.radians(45.0)   
#        h_weight = h
#        h_weight[ang_i > max_angle] = np.nan
#        
#        cor_elev = np.array([target.w_surf]).T - h_weight
#    
#        # append the "weighted" mean values for the actual depth and corrected elevation to
#        #   the target data frame   
#        target['h_avg_weight1'] = np.nanmean(h, axis = 1)
#        target['corElev_avg_weight1'] = np.nanmean(cor_elev, axis = 1)
#        
#        max_angle = np.radians(35.0)   
#        h_weight = h
#        h_weight[ang_i > max_angle] = np.nan
#        
#        cor_elev = np.array([target.w_surf]).T - h_weight
#    
#        # append the "weighted" mean values for the actual depth and corrected elevation to
#        #   the target data frame   
#        target['h_avg_weight2'] = np.nanmean(h, axis = 1)
#        target['corElev_avg_weight2'] = np.nanmean(cor_elev, axis = 1)
#        
#        max_angle = np.radians(30.0)   
#        h_weight = h
#        h_weight[ang_i > max_angle] = np.nan
#        
#        cor_elev = np.array([target.w_surf]).T - h_weight
#    
#        # append the "weighted" mean values for the actual depth and corrected elevation to
#        #   the target data frame   
#        target['h_avg_weight3'] = np.nanmean(h, axis = 1)
#        target['corElev_avg_weight3'] = np.nanmean(cor_elev, axis = 1)
    
    # return the target dataframe
#    return target, ang_r, x_dist, h, cor_elev
# END def correction
    
def pointFilter(tar,h,r,d):
    # tar = target out table from F(correction)
    # h = full corrected depth table from F(correction)
    # r = camera angles from F(visability)
    # d = camera distances from F(visability)
    
    max_ang = 35
    max_dist = 100
    
    h_filt = h
    h_filt[r > max_ang] = np.nan
    h_filt[d > max_dist] = np.nan
    
    tar['h_filt'] = np.nanmean(h_filt, axis = 1)
    
    cor_elev_filt = np.array([tar.w_surf]).T - h_filt
    tar['corElev_avg_filt'] = np.nanmean(cor_elev_filt, axis = 1)
    
    return tar

def timer(length,start_t):
    """timer function to calculate the running time"""
    
    num_proc = sum(length)    
    
    # time since processing started
    t_step = datetime.now() - start_t

    if t_step.total_seconds() <= 60:
        print("-> Finished %i points in %0.2f secs" %(num_proc,t_step.total_seconds()))
    else:
        ts = t_step.total_seconds() / 60    
        print("-> Finished %i points in %0.2f mins" %(num_proc,ts))
  
#END def timer
        
# END - MAIN PROGRAM HELPER FUNCTIONS

def main_prog():
    
    # CSV point cloud
    target_file =  'D:/Dropbox/Python/SfMAI/data/Mochlos_densesub_1-5mill.csv'
    
    # Camera Coords Exported from Metashape
    cam_file = 'D:/Dropbox/Python/SfMAI/data/mochlos_cams.csv'
    
    # sensor paramater file
    sensor_file = 'D:/Dropbox/Python/SfMAI/data/P3_sensor.csv'
    
    # output filename
    outfile = 'D:/Dropbox/Python/SfMAI/Mochlos_Test2.csv'
    
    # for a given dataset, the first run you can save the camera footprints (exportCam = True)
    #   for subsequent testing you can set exportCam = False, precalcCam = True
    #   pickle_file = path to exported pickle file
    exportCam = False
    precalcCam = True
    pickle_file = 'D:/Dropbox/Python/SfMAI/data/mochlos_cams_cam_foot.pkl' # .pkl
    

    extraOpt = np.array([True, False, False, True])

        
    print("The extra options are: ")
    print(extraOpt)
          
    # INPUTS - see sample_data folder in GitHub repository for file header formats
    
    # target points - as CSV point cloud (x,y,z,w_surf,r,g,b) from CloudCompare
    #   will be read in 10000 point chunks for memory managment purposes   
    targets = pd.read_csv(target_file, chunksize = 10000)

    
    # camera file - from Photoscan (Name, Position, Orientation...)
    # check for precalc checkbox, if so read directly to cam_r variable

    if precalcCam:
        foot_prints = pd.read_pickle(pickle_file)
        cams = pd.read_csv(cam_file)
    else:
        cams = pd.read_csv(cam_file)

    
    # camera sensor parameters - user generated csv
    sensor = pd.read_csv(sensor_file)
    
    # user feedback
    print("Data Loaded...")
    #sys.stdout.flush()
    
    # record the start time of the actual processing
    start_time = datetime.now()
    
    # array for count of total points
    count = []
    
    # Main Processing Loop, for each chunk of points from the reader
    for idx, tar in enumerate(targets):
        
        # Error check for mislabeled columns, from CloudCompare the column 
        #    header starts '//X,Y,Z
        if tar.columns.values[0] == '//X':
            tar.columns.values[0] = 'x'
            tar.columns.values[1] = 'y'
            tar.columns.values[2] = 'z'
        
        if tar.columns.values[0] == 'X':
            tar.columns.values[0] = 'x'
            tar.columns.values[1] = 'y'
            tar.columns.values[2] = 'z'
        
        count.append(tar.shape[0])    
        
        # if the index(idx) equals 0, use the mean elevation of the first chunk
        #   use the mean elevation of the first chunk of points to calculate the camera footprints
        if idx == 0 and not precalcCam:
            
            # establish mean elevation for footprint mapping from the mean 
            #   elevation of the target points
            base_elev = np.mean(tar.z)      
            
            # build camera footprints
            foot_prints = footprints(cams,sensor,base_elev)
            
            if exportCam:
                print("Saving footprints")
                foot_file = os.path.dirname(outfile) + '/' + os.path.basename(cam_file)[:-4]  + '_cam_foot.pkl'
                foot_prints.to_pickle(foot_file)
            
            # timer
            cam_end_time = datetime.now()
            mins_c = (cam_end_time - start_time).total_seconds() / 60
            print("Processed %i cameras in %0.2f minutes" %(np.count_nonzero(cams.x),mins_c))
        
        # use feedback and timer start
        if idx == 0 and precalcCam:
            refract_start_time = datetime.now()
        
        # test the visability of target point based on the camera footprints
        cam_r,cam_dist,cam_slpDist,radial_dist,cam_qual = visibility(cams,sensor,foot_prints,tar)
        
        # Save cam_R and for Debug targets
#        if self.exportCam_box.isChecked():
#        r_file = os.path.dirname(outfile) + '/' + os.path.basename(cam_file)[:-4]  + '_cam_r.csv'
#        np.savetxt(r_file,cam_r,delimiter=",")
#        tar_file = os.path.dirname(outfile) + '/d_tar.csv'
#        tar.to_csv(tar_file, header=True, index=False)
        
        
        # perform the refraction correction
#        tar_out, ang_r, x_dist, h, cor_elev = correction(cam_r, tar, extraOpt)
        
        # out dataframe
        tar_out = tar.copy(deep=True)
        
        #if filter checkbox is ticked
#        if extraOpt[3]:
#            tar_out = pointFilter(tar_out,h,cam_r,cam_dist)
            
        # OUTPUT STATS
        tar_out['cam_count'] = np.count_nonzero(~np.isnan(cam_r),axis=1)
            
        # QUALITY
        if "quality" in cams.columns:
            tar_out['cam_qual_mean'] = np.nanmean(cam_qual, axis = 1)
            tar_out['cam_qual_median'] = np.nanmean(cam_qual, axis = 1)
            tar_out['cam_qual_std'] = np.nanstd(cam_qual, axis = 1)
            tar_out['cam_qual_min'] = np.nanmin(cam_qual, axis = 1)
            tar_out['cam_qual_max'] = np.nanmax(cam_qual, axis = 1)
        else:
            tar_out['cam_qual_mean'] = 0
            tar_out['cam_qual_median'] = 0
            tar_out['cam_qual_std'] = 0
            tar_out['cam_qual_min'] = 0
            tar_out['cam_qual_max'] = 0
            
        # Angle
        tar_out['cam_ang_mean'] = np.nanmean(cam_r, axis = 1)
        tar_out['cam_ang_median'] = np.nanmean(cam_r, axis = 1)
        tar_out['cam_ang_std'] = np.nanstd(cam_r, axis = 1)
        tar_out['cam_ang_min'] = np.nanmin(cam_r, axis = 1)
        tar_out['cam_ang_max'] = np.nanmax(cam_r, axis = 1)
        
        # radial image dist
        tar_out['cam_rD_mean'] = np.nanmean(radial_dist, axis = 1)
        tar_out['cam_rD_median'] = np.nanmean(radial_dist, axis = 1)
        tar_out['cam_rD_std'] = np.nanstd(radial_dist, axis = 1)
        tar_out['cam_rD_min'] = np.nanmin(radial_dist, axis = 1)
        tar_out['cam_rD_max'] = np.nanmax(radial_dist, axis = 1)
        
        # Distance
        tar_out['cam_xydist_mean'] = np.nanmean(cam_dist, axis = 1)
        tar_out['cam_xydist_median'] = np.nanmean(cam_dist, axis = 1)
        tar_out['cam_xydist_std'] = np.nanstd(cam_dist, axis = 1)
        tar_out['cam_xydist_min'] = np.nanmin(cam_dist, axis = 1)
        tar_out['cam_xydist_max'] = np.nanmax(cam_dist, axis = 1)
        
        #slope distance
        tar_out['cam_slpDist_mean'] = np.nanmean(cam_slpDist, axis = 1)
        tar_out['cam_slpDist_median'] = np.nanmean(cam_slpDist, axis = 1)
        tar_out['cam_slpDist_std'] = np.nanstd(cam_slpDist, axis = 1)
        tar_out['cam_slpDist_min'] = np.nanmin(cam_slpDist, axis = 1)
        tar_out['cam_slpDist_max'] = np.nanmax(cam_slpDist, axis = 1)
        
        # DEBUG export
#        if extraOpt[2]:
#            file = os.path.dirname(outfile) + '/' + os.path.basename(target_file)[:-4]  + '_tar_out.pkl'
#            tar_out.to_pickle(file)
#            
#            file = os.path.dirname(outfile) + '/' + os.path.basename(target_file)[:-4]  + '_cam_r.csv'
#            np.savetxt(file,cam_r,delimiter=",")
#            
#            file = os.path.dirname(outfile) + '/' + os.path.basename(target_file)[:-4]  + '_cam_dist.csv'
#            np.savetxt(file,cam_dist,delimiter=",")
#            
#            file = os.path.dirname(outfile) + '/' + os.path.basename(target_file)[:-4]  + '_ang_r.csv'
#            np.savetxt(file,ang_r,delimiter=",")
#            
#            file = os.path.dirname(outfile) + '/' + os.path.basename(target_file)[:-4]  + '_x_dist.csv'
#            np.savetxt(file,x_dist,delimiter=",")
#            
#            file = os.path.dirname(outfile) + '/' + os.path.basename(target_file)[:-4]  + '_h_all.csv'
#            np.savetxt(file,h,delimiter=",")
#            
#            file = os.path.dirname(outfile) + '/' + os.path.basename(target_file)[:-4]  + '_cor_elev.csv'
#            np.savetxt(file,cor_elev,delimiter=",")
        
        # output - for the first chunk write header row, else append subsequent
        #   chunks without headers
        if idx == 0:
            with open(outfile, 'a') as f:
                tar_out.to_csv(f, header=True, index=False)
        else:
            with open(outfile, 'a') as f:
                 tar_out.to_csv(f, header=False, index=False)

        # user feedback, def timer and bottom progress bar
#        self.bot_progBar.setValue(idx)
        if precalcCam:
            timer(count, refract_start_time)
        else:
            timer(count, cam_end_time)
            
    # User feedback on the total processing time
    tot_count = sum(count)
    tot_time = (datetime.now() - start_time).total_seconds() / 60
    print("%i points processed, Total Running Time = %0.2f minutes" %(tot_count,tot_time))
#    self.botProg_lbl.setText('Processing Complete')

if __name__ == "__main__":
    main_prog()