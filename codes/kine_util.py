"""
(updated: 2017/03/02)

Class for kinematic analysis:
- plotting position / velocity / acceleration profiles
 - relative to body motion (or not)
 - smoothing via butterworth low-pass filter and spline curve fitting
- displaying animation for trajectories
 - or save as .png images and later combine to render .gif

"""
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import butter, filtfilt
from JSAnimation.IPython_display import display_animation, anim_to_html
from matplotlib import animation as animation
from scipy.interpolate import splrep, splev

class Kinematics(object):
    """
    Class computing and displaying the position/velocity/acceleration of motions of interest.
    """
    
    def __init__(self,vid_avo,vid_sho,avo_proj,sho_proj,
                 avo_pixl,sho_pixl,avo_fras,sho_fras,
                 avo_data,sho_data,cutoff):
        ## videos
        self.avo_vid = vid_avo
        self.sho_vid = vid_sho
        
        ## parameters
        self.avo_proj = avo_proj # projection angle
        self.sho_proj = sho_proj
        self.avo_pixl = avo_pixl # pixelwise physical length
        self.sho_pixl = sho_pixl
        self.avo_fras = avo_fras # framewise physical time
        self.sho_fras = sho_fras
        
        ## tracking data
        self.avo_data = avo_data
        self.sho_data = sho_data
        
        ### color coding for plotting uses
        self.clrs = ['r','b','g','y','c','m','w']
        
        ###
        self.cutoff = cutoff
        
    def disp_traj(self, showvid=True, option=1):
        """
        Animated trajectories overlayying (optional) on the videos.
        
         showvid == True will show the video; otherwise only points 
         
         option 1 = avocet
         option 2 = shoveler
        """
        
        ## which video
        if option == 1:
            vid = self.avo_vid
            xxxxs = [self.avo_data['pt1_cam1_X'], self.avo_data['pt2_cam1_X'],
                     self.avo_data['pt3_cam1_X'], self.avo_data['pt4_cam1_X'],
                     self.avo_data['pt5_cam1_X'], self.avo_data['pt6_cam1_X'], 
                     np.empty(vid.get_length())*np.nan]
            liney = [self.avo_data['pt1_cam1_Y'], self.avo_data['pt2_cam1_Y'],
                     self.avo_data['pt3_cam1_Y'], self.avo_data['pt4_cam1_Y'],
                     self.avo_data['pt5_cam1_Y'], self.avo_data['pt6_cam1_Y'], 
                     np.empty(vid.get_length())*np.nan]
        else:
            vid = self.sho_vid
            xxxxs = [self.sho_data['pt1_cam1_X'], self.sho_data['pt2_cam1_X'],
                     self.sho_data['pt3_cam1_X'], self.sho_data['pt4_cam1_X'],
                     self.sho_data['pt5_cam1_X'], self.sho_data['pt6_cam1_X'], 
                     self.sho_data['pt7_cam1_X']]
            liney = [self.sho_data['pt1_cam1_Y'], self.sho_data['pt2_cam1_Y'],
                     self.sho_data['pt3_cam1_Y'], self.sho_data['pt4_cam1_Y'],
                     self.sho_data['pt5_cam1_Y'], self.sho_data['pt6_cam1_Y'], 
                     self.sho_data['pt7_cam1_Y']]
        clrs = self.clrs
        
        ##
        totalframes = vid.get_length()
        
        fig = plt.figure(figsize=(24,16))
        ax  = plt.axes(xlim=(0, vid.get_data(0).shape[1]), ylim=(vid.get_data(0).shape[0],0))
        
        lines = sum([ax.plot([], [], '-', c=c, lw=10) for c in clrs], [])
        
        def init():
            for line in lines:
                line.set_data([], [])
            return lines
        
        def animate(i):
            # video
            if showvid:
                image = vid.get_data(i)
                ax.imshow(image)
            
            j = 0
            for line in lines:
                line.set_data(xxxxs[j][:i],liney[j][:i])
                j += 1
            return lines
        ani = animation.FuncAnimation(fig, animate, init_func=init,
                                      frames=totalframes, repeat=False)
        return ani
        
    def traj_png(self, option=1):
        """
        Since the "disp_traj" method is too slow, 
        each frame is saved as a png image in this method.
        
        Additional features (compared to "disp_traj") are
        - figure title
        - 
        """
        
        ## which video
        if option == 1:
            vid   = self.avo_vid
            time  = self.avo_fras
            name  = 'Avo'
            xxxxs = [self.avo_data['pt1_cam1_X'], self.avo_data['pt2_cam1_X'],
                     self.avo_data['pt3_cam1_X'], self.avo_data['pt4_cam1_X'],
                     self.avo_data['pt5_cam1_X'], self.avo_data['pt6_cam1_X'], 
                     np.empty(vid.get_length())*np.nan]
            liney = [self.avo_data['pt1_cam1_Y'], self.avo_data['pt2_cam1_Y'],
                     self.avo_data['pt3_cam1_Y'], self.avo_data['pt4_cam1_Y'],
                     self.avo_data['pt5_cam1_Y'], self.avo_data['pt6_cam1_Y'], 
                     np.empty(vid.get_length())*np.nan]
        else:
            vid   = self.sho_vid
            time  = self.sho_fras
            name  = 'Sho'
            xxxxs = [self.sho_data['pt1_cam1_X'], self.sho_data['pt2_cam1_X'],
                     self.sho_data['pt3_cam1_X'], self.sho_data['pt4_cam1_X'],
                     self.sho_data['pt5_cam1_X'], self.sho_data['pt6_cam1_X'], 
                     self.sho_data['pt7_cam1_X']]
            liney = [self.sho_data['pt1_cam1_Y'], self.sho_data['pt2_cam1_Y'],
                     self.sho_data['pt3_cam1_Y'], self.sho_data['pt4_cam1_Y'],
                     self.sho_data['pt5_cam1_Y'], self.sho_data['pt6_cam1_Y'], 
                     self.sho_data['pt7_cam1_Y']]
        clrs  = self.clrs
        label = ['RW','Head','Toe','Neck (base)','Tail (base)','LW','Ref.']
        
        ## dimensions
        xlim = vid.get_data(0).shape[1]
        ylim = vid.get_data(0).shape[0]
        
        ## loop through frames
        for frame in range(vid.get_length()):
            plt.figure(figsize=(24,16))
            
            ## video frame
            image = vid.get_data(frame)
            plt.imshow(image)
            
            ## trajectories
            for j in range(len(liney)):
                plt.plot(xxxxs[j][:frame],liney[j][:frame],c=clrs[j],lw=8,label=label[j])
                
            plt.title('Time: %.3f s' %(time*frame), fontsize=32)
            plt.xlim(0,xlim); plt.ylim(ylim,0)
            plt.legend(loc='upper right', fontsize=24)

            plt.savefig('./%s_traj/frame_%03d.png' %(name,frame))
            plt.close() # prevent displaying inline
        
    def vel_acc(self):
        """
        Compute the velocity and acceleration profiles.
        """
        ### read in and correct for 'global' motions (camera)
        ## camera
        sho_cax,sho_cay = np.array(self.sho_data['pt7_cam1_X']), \
                          np.array(self.sho_data['pt7_cam1_Y'])
        
        ## wing tips
        avo_rwx,avo_rwy = np.array(self.avo_data['pt1_cam1_X']), \
                          np.array(self.avo_data['pt1_cam1_Y'])
        sho_rwx,sho_rwy = np.array(self.sho_data['pt1_cam1_X'])-sho_cax, \
                         -np.array(self.sho_data['pt1_cam1_Y'])+sho_cay
        avo_lwx,avo_lwy = np.array(self.avo_data['pt6_cam1_X']), \
                          np.array(self.avo_data['pt6_cam1_Y'])
        sho_lwx,sho_lwy = np.array(self.sho_data['pt6_cam1_X'])-sho_cax, \
                         -np.array(self.sho_data['pt6_cam1_Y'])+sho_cay
        
        ## head (beak or eye)
        avo_hdx,avo_hdy = np.array(self.avo_data['pt2_cam1_X']), \
                          np.array(self.avo_data['pt2_cam1_Y'])
        sho_hdx,sho_hdy = np.array(self.sho_data['pt2_cam1_X'])-sho_cax, \
                         -np.array(self.sho_data['pt2_cam1_Y'])+sho_cay
        
        ## toe
        avo_tox,avo_toy = np.array(self.avo_data['pt3_cam1_X']), \
                          np.array(self.avo_data['pt3_cam1_Y'])
        sho_tox,sho_toy = np.array(self.sho_data['pt3_cam1_X'])-sho_cax, \
                         -np.array(self.sho_data['pt3_cam1_Y'])+sho_cay
        
        ## body (neck and tail bases)
        avo_nex,avo_ney = np.array(self.avo_data['pt4_cam1_X']), \
                          np.array(self.avo_data['pt4_cam1_Y'])
        sho_nex,sho_ney = np.array(self.sho_data['pt4_cam1_X'])-sho_cax, \
                         -np.array(self.sho_data['pt4_cam1_Y'])+sho_cay
        avo_tax,avo_tay = np.array(self.avo_data['pt5_cam1_X']), \
                          np.array(self.avo_data['pt5_cam1_Y'])
        sho_tax,sho_tay = np.array(self.sho_data['pt5_cam1_X'])-sho_cax, \
                         -np.array(self.sho_data['pt5_cam1_Y'])+sho_cay
        
        ##- body motion correction
        avo_rwcx,avo_rwcy = avo_rwx-avo_nex,avo_ney-avo_rwy
        avo_hdcx,avo_hdcy = avo_hdx-avo_nex,avo_ney-avo_hdy
        avo_tocx,avo_tocy = avo_tox-avo_nex,avo_ney-avo_toy
        avo_necx,avo_necy = avo_nex-avo_nex,avo_ney-avo_ney
        avo_tacx,avo_tacy = avo_tax-avo_nex,avo_ney-avo_tay
        avo_lwcx,avo_lwcy = avo_lwx-avo_nex,avo_ney-avo_lwy
        
        sho_rwcx,sho_rwcy = sho_rwx-sho_nex,sho_ney-sho_rwy
        sho_hdcx,sho_hdcy = sho_hdx-sho_nex,sho_ney-sho_hdy
        sho_tocx,sho_tocy = sho_tox-sho_nex,sho_ney-sho_toy
        sho_necx,sho_necy = sho_nex-sho_nex,sho_ney-sho_ney
        sho_tacx,sho_tacy = sho_tax-sho_nex,sho_ney-sho_tay
        sho_lwcx,sho_lwcy = sho_lwx-sho_nex,sho_ney-sho_lwy

        
        ## wrap up and calibrate for physical sizes
        avo_x  = [avo_rwx*self.avo_pixl ,avo_hdx*self.avo_pixl ,avo_tox*self.avo_pixl, \
                  avo_nex*self.avo_pixl ,avo_tax*self.avo_pixl ,avo_lwx*self.avo_pixl]
        avo_cx = [avo_rwcx*self.avo_pixl,avo_hdcx*self.avo_pixl,avo_tocx*self.avo_pixl, \
                  avo_necx*self.avo_pixl,avo_tacx*self.avo_pixl,avo_lwcx*self.avo_pixl]
        avo_y  = [avo_rwy*self.avo_pixl ,avo_hdy*self.avo_pixl ,avo_toy*self.avo_pixl, \
                  avo_ney*self.avo_pixl ,avo_tay*self.avo_pixl ,avo_lwy*self.avo_pixl]
        avo_cy = [avo_rwcy*self.avo_pixl,avo_hdcy*self.avo_pixl,avo_tocy*self.avo_pixl, \
                  avo_necy*self.avo_pixl,avo_tacy*self.avo_pixl,avo_lwcy*self.avo_pixl]
        sho_x  = [sho_rwx*self.sho_pixl ,sho_hdx*self.sho_pixl ,sho_tox*self.sho_pixl, \
                  sho_nex*self.sho_pixl ,sho_tax*self.sho_pixl ,sho_lwx*self.sho_pixl]
        sho_cx = [sho_rwcx*self.sho_pixl,sho_hdcx*self.sho_pixl,sho_tocx*self.sho_pixl, \
                  sho_necx*self.sho_pixl,sho_tacx*self.sho_pixl,sho_lwcx*self.sho_pixl]
        sho_y  = [sho_rwy*self.sho_pixl ,sho_hdy*self.sho_pixl ,sho_toy*self.sho_pixl, \
                  sho_ney*self.sho_pixl ,sho_tay*self.sho_pixl ,sho_lwy*self.sho_pixl]
        sho_cy = [sho_rwcy*self.sho_pixl,sho_hdcy*self.sho_pixl,sho_tocy*self.sho_pixl, \
                  sho_necy*self.sho_pixl,sho_tacy*self.sho_pixl,sho_lwcy*self.sho_pixl]
        
        ### time stamp correction
        avo_time = np.linspace(0.,len(avo_rwx)*self.avo_fras,len(avo_rwx))
        sho_time = np.linspace(0.,len(sho_rwx)*self.sho_fras,len(sho_rwx))
        
        ### submethods to fit spline curves
        def _spl_fit(t,x,y,species,c=False,k=3,sl=5,sm=5):
            rwx,hdx,tox,nex,tax,lwx = x
            rwy,hdy,toy,ney,tay,lwy = y
            
            ## "un-nan" the data
            ## - rotation angle
            if c is True:
                where = np.argwhere(~np.isnan(nex)&~np.isnan(tax))
                ro_t = t[where][:,0]
                ro_p = np.rad2deg(np.arctan(tay/tax))[where][:,0]
            else:
                ro_t,ro_p=np.zeros(100),np.zeros(100)
            
            ## - times and original positions
            tt = [self._unnan(t,rwx)[0],self._unnan(t,hdx)[0],self._unnan(t,tox)[0], \
                  self._unnan(t,nex)[0],self._unnan(t,tax)[0],self._unnan(t,lwx)[0]]
            xx = [self._unnan(t,rwx)[1],self._unnan(t,hdx)[1],self._unnan(t,tox)[1], \
                  self._unnan(t,nex)[1],self._unnan(t,tax)[1],self._unnan(t,lwx)[1]]
            yy = [self._unnan(t,rwy)[1],self._unnan(t,hdy)[1],self._unnan(t,toy)[1], \
                  self._unnan(t,ney)[1],self._unnan(t,tay)[1],self._unnan(t,lwy)[1]]
            
            
            ## smoothing
            x_,y_ = [],[]
            for i in range(len(tt)):
                x_.append(self._low_pass(data=xx[i],species=species))
                y_.append(self._low_pass(data=yy[i],species=species))
            ro_p = self._low_pass(ro_p,species=species)
            
            ## fit splines
            rwx_f,rwy_f = splrep(tt[0],x_[0],k=k,s=sl), splrep(tt[0],y_[0],k=k,s=sl)
            hdx_f,hdy_f = splrep(tt[1],x_[1],k=k,s=sm), splrep(tt[1],y_[1],k=k,s=sm)
            tox_f,toy_f = splrep(tt[2],x_[2],k=k,s=sl), splrep(tt[2],y_[2],k=k,s=sl)
            nex_f,ney_f = splrep(tt[3],x_[3],k=k,s=sl), splrep(tt[3],y_[3],k=k,s=sl)
            tax_f,tay_f = splrep(tt[4],x_[4],k=k,s=sm), splrep(tt[4],y_[4],k=k,s=sm)
            lwx_f,lwy_f = splrep(tt[5],x_[5],k=k,s=sl), splrep(tt[5],y_[5],k=k,s=sl)
            ro_f        = splrep(ro_t ,ro_p ,k=k,s=sl)
        
                        
            ## - positions
            rwpx_f,rwpy_f = splev(tt[0],rwx_f),splev(tt[0],rwy_f)
            hdpx_f,hdpy_f = splev(tt[1],hdx_f),splev(tt[1],hdy_f)
            topx_f,topy_f = splev(tt[2],tox_f),splev(tt[2],toy_f)
            nepx_f,nepy_f = splev(tt[3],nex_f),splev(tt[3],ney_f)
            tapx_f,tapy_f = splev(tt[4],tax_f),splev(tt[4],tay_f)
            lwpx_f,lwpy_f = splev(tt[5],lwx_f),splev(tt[5],lwy_f)
            ro_pf         = splev(ro_t ,ro_f)
            
            px = [rwpx_f,hdpx_f,topx_f,nepx_f,tapx_f,lwpx_f]
            py = [rwpy_f,hdpy_f,topy_f,nepy_f,tapy_f,lwpy_f]
            
            ## - velocities
            rwvx_f,rwvy_f = splev(tt[0],rwx_f,der=1),splev(tt[0],rwy_f,der=1)
            hdvx_f,hdvy_f = splev(tt[1],hdx_f,der=1),splev(tt[1],hdy_f,der=1)
            tovx_f,tovy_f = splev(tt[2],tox_f,der=1),splev(tt[2],toy_f,der=1)
            nevx_f,nevy_f = splev(tt[3],nex_f,der=1),splev(tt[3],ney_f,der=1)
            tavx_f,tavy_f = splev(tt[4],tax_f,der=1),splev(tt[4],tay_f,der=1)
            lwvx_f,lwvy_f = splev(tt[5],lwx_f,der=1),splev(tt[5],lwy_f,der=1)
            ro_vf         = splev(ro_t ,ro_f ,der=1)
            
            vx = [rwvx_f,hdvx_f,tovx_f,nevx_f,tavx_f,lwvx_f]
            vy = [rwvy_f,hdvy_f,tovy_f,nevy_f,tavy_f,lwvy_f]
            
            ## - accelerations
            rwax_f,rway_f = splev(tt[0],rwx_f,der=2),splev(tt[0],rwy_f,der=2)
            hdax_f,hday_f = splev(tt[1],hdx_f,der=2),splev(tt[1],hdy_f,der=2)
            toax_f,toay_f = splev(tt[2],tox_f,der=2),splev(tt[2],toy_f,der=2)
            neax_f,neay_f = splev(tt[3],nex_f,der=2),splev(tt[3],ney_f,der=2)
            taax_f,taay_f = splev(tt[4],tax_f,der=2),splev(tt[4],tay_f,der=2)
            lwax_f,lway_f = splev(tt[5],lwx_f,der=2),splev(tt[5],lwy_f,der=2)
            ro_af         = splev(ro_t ,ro_f ,der=2)
            
            ax = [rwax_f,hdax_f,toax_f,neax_f,taax_f,lwax_f]
            ay = [rway_f,hday_f,toay_f,neay_f,taay_f,lway_f]

            ro = [ro_t,ro_pf,ro_vf,ro_af]
            
            return tt,xx,yy,px,py,vx,vy,ax,ay,ro
        
        k, sl, sm = 5, 3, 3
        avo_tpva = _spl_fit(avo_time,avo_x,avo_y,species="avo")
        sho_tpva = _spl_fit(sho_time,sho_x,sho_y,species="sho")
        avo_cpva = _spl_fit(avo_time,avo_cx,avo_cy,species="avo",c=True)
        sho_cpva = _spl_fit(sho_time,sho_cx,sho_cy,species="sho",c=True)
        
        return avo_tpva,sho_tpva,avo_cpva,sho_cpva
        
    def disp_vel_acc(self):
        """
        Display the velocity and acceleration profiles
        """
        ### extract original/fitted data
        avo_tpva, sho_tpva, avo_cpva, sho_cpva = self.vel_acc()
                
        ### plotting
        def _plot(tpva,n):
            t,x,y,px,py,vx,vy,ax,ay,ro = tpva
            rot,rop,rov,roa = ro
            
            plt.figure(figsize=(24,70))
            ## frame-moving/body-motion (def:neck base) corrected positions
            plt.subplot(421)
            plt.plot(t[0],x[0],'-.',c=self.clrs[0],lw=5); plt.plot(t[0],px[0],c=self.clrs[0],label='RW')
            plt.plot(t[1],x[1],'-.',c=self.clrs[1],lw=5); plt.plot(t[1],px[1],c=self.clrs[1],label='Head')
            plt.plot(t[2],x[2],'-.',c=self.clrs[2],lw=5); plt.plot(t[2],px[2],c=self.clrs[2],label='Toe')
            plt.plot(t[3],x[3],'-.',c=self.clrs[3],lw=5); plt.plot(t[3],px[3],c=self.clrs[3],label='Neck (base)')
            plt.plot(t[4],x[4],'-.',c=self.clrs[4],lw=5); plt.plot(t[4],px[4],c=self.clrs[4],label='Tail (base)')
            plt.plot(t[5],x[5],'-.',c=self.clrs[5],lw=5); plt.plot(t[5],px[5],c=self.clrs[5],label='LW')
            plt.title('X Positions [%s]' %n, fontsize=32); plt.legend(fontsize=20, loc='upper left')
            plt.xlabel('Time (s)',fontsize=24); plt.ylabel('Position (cm)',fontsize=24)
            plt.xticks(fontsize=20), plt.yticks(fontsize=20)
            
            plt.subplot(422)
            plt.plot(t[0],y[0],'-.',c=self.clrs[0],lw=5); plt.plot(t[0],py[0],c=self.clrs[0],label='RW')
            plt.plot(t[1],y[1],'-.',c=self.clrs[1],lw=5); plt.plot(t[1],py[1],c=self.clrs[1],label='Head')
            plt.plot(t[2],y[2],'-.',c=self.clrs[2],lw=5); plt.plot(t[2],py[2],c=self.clrs[2],label='Toe')
            plt.plot(t[3],y[3],'-.',c=self.clrs[3],lw=5); plt.plot(t[3],py[3],c=self.clrs[3],label='Neck (base)')
            plt.plot(t[4],y[4],'-.',c=self.clrs[4],lw=5); plt.plot(t[4],py[4],c=self.clrs[4],label='Tail (base)')
            plt.plot(t[5],y[5],'-.',c=self.clrs[5],lw=5); plt.plot(t[5],py[5],c=self.clrs[5],label='LW')
            plt.title('Y Positions [%s]' %n, fontsize=32); plt.legend(fontsize=20, loc='upper left')
            plt.xlabel('Time (s)',fontsize=24); plt.ylabel('Position (cm)',fontsize=24)
            plt.xticks(fontsize=20), plt.yticks(fontsize=20)
            
            ## velocities
            plt.subplot(423)
            plt.plot(t[0],vx[0],c=self.clrs[0],label='RW')
            plt.plot(t[1],vx[1],c=self.clrs[1],label='Head')
            plt.plot(t[2],vx[2],c=self.clrs[2],label='Toe')
            plt.plot(t[3],vx[3],c=self.clrs[3],label='Neck (base)')
            plt.plot(t[4],vx[4],c=self.clrs[4],label='Tail (base)')
            plt.plot(t[5],vx[5],c=self.clrs[5],label='LW')
            plt.title('X Velocities [%s]' %n, fontsize=32); plt.legend(fontsize=20, loc='upper right')
            plt.xlabel('Time (s)',fontsize=24); plt.ylabel('Velocity (cm/s)',fontsize=24)
            plt.xticks(fontsize=20), plt.yticks(fontsize=20); plt.ylim(-1500,1500)
            
            plt.subplot(424)
            plt.plot(t[0],vy[0],c=self.clrs[0],label='RW')
            plt.plot(t[1],vy[1],c=self.clrs[1],label='Head')
            plt.plot(t[2],vy[2],c=self.clrs[2],label='Toe')
            plt.plot(t[3],vy[3],c=self.clrs[3],label='Neck (base)')
            plt.plot(t[4],vy[4],c=self.clrs[4],label='Tail (base)')
            plt.plot(t[5],vy[5],c=self.clrs[5],label='LW')
            plt.title('Y Velocities [%s]' %n, fontsize=32); plt.legend(fontsize=20, loc='upper right')
            plt.xlabel('Time (s)',fontsize=24); plt.ylabel('Velocity (cm/s)',fontsize=24)
            plt.xticks(fontsize=20), plt.yticks(fontsize=20); plt.ylim(-1500,1500)
            
            ## accelerations
            plt.subplot(425)
            plt.plot(t[0],ax[0],c=self.clrs[0],label='RW')
            plt.plot(t[1],ax[1],c=self.clrs[1],label='Head')
            plt.plot(t[2],ax[2],c=self.clrs[2],label='Toe')
            plt.plot(t[3],ax[3],c=self.clrs[3],label='Neck (base)')
            plt.plot(t[4],ax[4],c=self.clrs[4],label='Tail (base)')
            plt.plot(t[5],ax[5],c=self.clrs[5],label='LW')
            plt.title('X Accelerations [%s]' %n, fontsize=32); plt.legend(fontsize=20, loc='upper right')
            plt.xlabel('Time (s)',fontsize=24); plt.ylabel('Acceleration (cm/s^2)',fontsize=24)
            plt.xticks(fontsize=20), plt.yticks(fontsize=20); plt.ylim(-15000,15000)
            
            plt.subplot(426)
            plt.plot(t[0],ay[0],c=self.clrs[0],label='RW')
            plt.plot(t[1],ay[1],c=self.clrs[1],label='Head')
            plt.plot(t[2],ay[2],c=self.clrs[2],label='Toe')
            plt.plot(t[3],ay[3],c=self.clrs[3],label='Neck (base)')
            plt.plot(t[4],ay[4],c=self.clrs[4],label='Tail (base)')
            plt.plot(t[5],ay[5],c=self.clrs[5],label='LW')
            plt.title('Y Accelerations [%s]' %n, fontsize=32); plt.legend(fontsize=20, loc='upper right')
            plt.xlabel('Time (s)',fontsize=24); plt.ylabel('Acceleration (cm/s^2)',fontsize=24)
            plt.xticks(fontsize=20), plt.yticks(fontsize=20); plt.ylim(-15000,15000)
            
            ## rotation
            ## = arcTan[( neck - tail )y / ( neck - tail )x ]
            plt.subplot(427)
            plt.plot(rot,rop,c=self.clrs[0],label='Angle')
            plt.plot(rot,rov/10,c=self.clrs[1],label='Angular Vel. (/10)')
            plt.plot(rot,roa/1000,c=self.clrs[2],label='Angular Acc. (/1000)')
            plt.title('Body Rotation (0=horizontal) [%s]' %n, fontsize=32); plt.legend(fontsize=20, loc='upper right')
            plt.xlabel('Time (s)',fontsize=24); plt.ylabel('Angle (deg) [velocity/acceleration]',fontsize=24)
            plt.xticks(fontsize=20), plt.yticks(fontsize=20); plt.ylim(-60,60)
        
        _plot(avo_tpva,n='Avocet')
        _plot(sho_tpva,n='Shoveler')
        _plot(avo_cpva,n='Avocet [rel.]')
        _plot(sho_cpva,n='Shoveler [rel.]')
        
    def coord_trans(self):
        """
        Transforming 'inertial' frame (camera) to 'rotating' frame (bird)
        """
        pass
    
    #########################
    # Supplementary methods #
    def _unnan(self,t,data):
        """
        Return the non-nan sections of the data and the time stamps
        """
        ix = np.argwhere(~np.isnan(data))
        t_ = t[ix]
        d_ = data[ix]
        
        return t_[:,0], d_[:,0]
    
    
    def _low_pass(self,data,species,order=5):
        """
        Butterworth low-pass filter
        
        Parameters
         fs: float
             Sampling frequency. Default to 1 / frame rate
         cutoff: float
             Cutoff level as a multiple of fs / 100.
             Estimated from the curlew video [youtu.be/Zc6P8LaiJco]
             which shows ~5Hz flapping frequency, therefore 5 = 500 / 100.
             (Note: Curlew has similar wing span)
             Default to 3 times of the flapping frequency.
         order: integer
             The order of the Butterworth filter
        """
        ##
        cutoff = self.cutoff
        if species=="avo":
            fs = self.avo_fras
        else:
            fs = self.sho_fras
        
        ##
        def __butter_lowpass(cutoff, fs, order=order):
            nyq = 0.5 * fs
            normal_cutoff = cutoff / nyq
            b, a = butter(order, normal_cutoff, btype='low', analog=False)
            return b, a

        def __butter_lowpass_filter(data, cutoff, fs, order=order):
            b, a = __butter_lowpass(cutoff, fs, order)
            y = filtfilt(b, a, data)
            return y
        
        ## shift to zero-starting
        ini_d = data[0]
        data_ = data - ini_d
        
        ## 
        data_filt = __butter_lowpass_filter(data=data_,cutoff=fs/100*cutoff,fs=fs)
        
        ## restore the translation
        return data_filt + ini_d
