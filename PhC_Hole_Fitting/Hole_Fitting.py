import cv2
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
import scipy.io as sio
import numpy as np
import math
import os
from operator import add


def gauss(x, amp, sigma, x0):
        return amp * np.exp(-(x-x0)**2/(2*sigma**2))

fol = "dev" #folder containing all the SEM images of a given device. Here, "dev" = "device"
dev_num = [6] #The devices to analyze
part_num = np.linspace(1,13,13,dtype=np.int) #The number of images within a given device folder
nm_per_px = 200/518 #nanometres per pixel. This can be found from the scale bar in the SEM image

for m in range(len(dev_num)):
    for i in range(len(part_num)-1):
        fname = "holes_0"
        if i < 9:
            fname = "holes_00"
        fpath = fol + str(dev_num[m]) + "/" + fname + str(part_num[i]) + ".jpg"
        img = cv2.imread(fpath)

        # Convert into grayscale and threshed it
        threshold = 100

        gray = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)
        th, threshed = cv2.threshold(gray, threshold, 255, cv2.THRESH_BINARY)

        # Morph to denoise
        threshed = cv2.dilate(threshed, None)
        threshed = cv2.erode(threshed, None)

        # Find the external contours
        cnts,hierarchy = cv2.findContours(threshed, cv2.RETR_CCOMP, cv2.CHAIN_APPROX_NONE)
        inner_cnts = [cnts[l] for l in range(len(cnts)) if hierarchy[0,l,3] >= 0]
        holes = [inner_cnts[l] for l in range(len(inner_cnts)) if len(inner_cnts[l]) > 200 
                 and cv2.contourArea(inner_cnts[l]) > 30000
                 and len(inner_cnts[l]) < 4000]

        c_img = cv2.drawContours(img, holes, -1, (255, 0, 0), 2, cv2.LINE_AA)

        # Fit ellipses
        fitResult = []
        for cnt in holes:
            rbox = cv2.fitEllipse(cnt)
            fitResult.append(rbox)
            cv2.ellipse(img, rbox, (255, 100, 255), 2, cv2.LINE_AA)

        cv2.imwrite(fol + str(dev_num[m]) + "/" + fname + str(dev_num[m]) + "_" + str(part_num[i]) + "to" + str(part_num[i+1]) + "_fitted.jpg",img)
        img2 = img[:,:,::-1]

        hole_xcor = []
        hole_ycor = []
        hole_hx = []
        hole_hy = []

        hole_xcor = [nm_per_px*fitResult[l][0][0] for l in range(len(fitResult))]
        hole_ycor = [nm_per_px*fitResult[l][0][1] for l in range(len(fitResult))]
        hole_hx = [nm_per_px*fitResult[l][1][0] for l in range(len(fitResult))]
        hole_hy = [nm_per_px*fitResult[l][1][1] for l in range(len(fitResult))]

        print([len(holes),len(fitResult)])

        if len(holes) == len(fitResult):
            angle = [fitResult[l][2]/180*math.pi for l in range(len(fitResult))]

            sigma_list = []
            x0_list = []
            fig, axs = plt.subplots(len(fitResult))
            for k in range(len(fitResult)):
                diff_nm_list = []
                for j in range(len(holes[k])):
                    posx = holes[k][j][0][0] - fitResult[k][0][0]
                    posy = holes[k][j][0][1]- fitResult[k][0][1]
                    r = math.sqrt(posx*posx+posy*posy)
                    theta = math.atan(posy/posx)
                    r_expected = 0.5*fitResult[k][1][0]*fitResult[k][1][1]/math.sqrt(
                        fitResult[k][1][1]*math.cos(theta)*fitResult[k][1][1]*math.cos(theta)
                        +fitResult[k][1][0]*math.sin(theta)*fitResult[k][1][0]*math.sin(theta))
                    diff_nm = (r_expected - r)*nm_per_px
                    diff_nm_list.append(diff_nm)

                data_entries, bins = np.histogram(diff_nm_list,bins=np.linspace(-10, 10, 401))
                binscenters = np.array([0.5 * (bins[l] + bins[l+1]) for l in range(len(bins)-1)])
                popt, pcov = curve_fit(gauss, xdata=binscenters, ydata=data_entries)

                if len(fitResult) == 1:
                    axs.hist(diff_nm_list,bins=np.linspace(-10, 10, 401))
                    axs.plot(binscenters,gauss(binscenters,*popt))
                else:
                    axs[k].hist(diff_nm_list,bins=np.linspace(-10, 10, 401))
                    axs[k].plot(binscenters,gauss(binscenters,*popt))

                x0_list.append(popt[2])
                sigma_list.append(popt[1])

            fig.savefig(fol + str(dev_num[m]) + "/" + fname + str(dev_num[m]) + "_" + str(part_num[i]) + "to" + str(part_num[i+1]) + "_FitResult.jpg")
            fig.clear()
            plt.close(fig)
            print(sigma_list)
            print(x0_list)
            hole_hx_modified = list(map(add, hole_hx, x0_list))
            hole_hy_modified = list(map(add, hole_hy, x0_list))

            if os.path.isfile(fol + str(dev_num[m]) + "/" + fname+ str(part_num[i]) + "to" + str(part_num[i+1]) +  '.mat'):
                os.remove(fol + str(dev_num[m]) + "/" + fname + str(part_num[i]) + "to" + str(part_num[i+1]) +  '.mat')
            sio.savemat(fol + str(dev_num[m]) + "/" + str(part_num[i]) + "to" + str(part_num[i+1]) +  '.mat',
                        {'sigma':sigma_list, 
                         'hole_xcor':hole_xcor, 
                         'hole_ycor':hole_ycor, 
                         'hole_hx':hole_hx_modified, 
                         'hole_hy':hole_hy_modified,
                         'angle':angle})