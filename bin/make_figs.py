#!/usr/bin/env python
# encoding: utf-8
"""
untitled.py

Created by labuser on 2010-10-20.
Copyright (c) 2010 __MyCompanyName__. All rights reserved.
"""

import sys
import os
import cPickle as pk
import pylab as pl
import numpy as np

def adjust_spines(ax,spines):
    for loc, spine in ax.spines.iteritems():
        if loc in spines:
            spine.set_position(('outward',10)) # outward by 10 points
            #spine.set_smart_bounds(True)
        else:
            spine.set_color('none') # don't draw spine
    # turn off ticks where there is no spine
    if 'left' in spines:
        ax.yaxis.set_ticks_position('left')
    else:
        # no yaxis ticks
        ax.yaxis.set_ticks([])
    if 'bottom' in spines:
        ax.xaxis.set_ticks_position('bottom')
    else:
        # no xaxis ticks
        ax.xaxis.set_ticks([])


def get_fp(data):
	stim_ind_all = data[17]['b_sp_im_ind']
	stim_id_all = data[17]['b_sp_im_id']
	stim_ind = np.unique(np.array(data[17]['b_sp_im_ind']))
	fp = {}
	for ch in [3,17,27,37,41,42,43,45,48,49,51,52,53,55,56,58,69,70,75,79,83,84,85,89]:#data.keys()[:50]:
		if data[ch]['good_sp'] != []:
			if ch not in fp.keys():
				fp[ch] = {}
			ch_stim_ind = np.array(data[ch]['g_sp_im_ind'])
			sp_t = np.array(data[ch]['g_sp_t_rel'])
			for s_ind in stim_ind:
				t_im_id = int(stim_id_all[stim_ind_all.index(s_ind)].split('_')[1])
				if t_im_id not in fp[ch].keys():
					fp[ch][t_im_id] = []
				sp_in = np.flatnonzero(ch_stim_ind == s_ind)
				if len(sp_in) == 0:
					fp[ch][t_im_id].append(0)
				else:
					c_sp_t = sp_t[sp_in]
					fp[ch][t_im_id].append(len(c_sp_t > 75000))
	fp_fr = np.empty((len(fp),len(fp[17].keys())))
	for ch_num,ch in enumerate(fp.keys()):
		for im_num,im_ind in enumerate(fp[ch].keys()):
			fp_fr[ch_num,im_num] = np.mean(np.array(fp[ch][im_ind]))/0.175 
	for p_num in range(24):
		ax = pl.subplot(5,5,p_num+1)
		ax.bar(range(1,106),fp_fr[p_num,:].T,color = [0.5,0.5,0.5],edgecolor=  [0.5,0.5,0.5])
		ax.set_xlim([1,105])
		#print max(fp_fr)
		ax.set_ylim([0,np.max(fp_fr)])
		if p_num == 20:
			adjust_spines(ax,['left','bottom'])
			ax.set_xlabel('Image #')
			ax.set_ylabel('Response (ips)')
			ax.set_xticks([0,50,100])
			ax.set_xticklabels(['0','50','100'])
		else:
			adjust_spines(ax,[])
	pl.savefig('Bar_plots_24.pdf')		
			
			
			
			 
def cluster_neurons(data):
	fp = []
	ch_num = []
	for ch in data.keys():
		if data[ch]['good_sp'] != []:
			stim_id = [int(foo.split('_')[1]) for foo in data[ch]['g_sp_im_id']]
			n,b = np.histogram(stim_id,bins=range(0,106,1))
			fp.append(n[:100])
			ch_num.append(ch)
	ind_seed = 38
	fp_c = fp[ind_seed]
	ch_o = []
	ch_o.append(ch_num[ind_seed])
	fp[ind_seed] = np.zeros(np.shape(fp[ind_seed]))
	while sum(sum(fp)) > 0:
		n_fp = np.argmax(np.dot(np.array(fp),fp_c))
		ch_o.append(ch_num[n_fp])
		fp_c = fp[n_fp]
		fp[n_fp] = np.zeros(np.shape(fp[n_fp]))
	return ch_o
		 
		
		
def plot_crude_fp(data,channels):
	sp_ind = 1
	for ch in channels:
		pl.subplot(10,4,sp_ind)
		if data[ch]['good_sp'] != []:
			pl.hist([int(foo.split('_')[1]) for foo in data[ch]['g_sp_im_id']],bins = range(1,105,1))
			sp_ind +=1
		if sp_ind > 40:
			break
	pl.show()
	
def plot_raster(data,imdir):
	for first_img_ind in [50]:#[1,10,30,40,50,100,150,200,250,300,400,500,600,650,800,1000]:
		pl.figure()
		pl.clf()
		for c_im_ind in range(1,6,1): 
			ch_num = 1
			ax = pl.subplot(2,5,c_im_ind+5)
			im = []
			for ch in data.keys():
				if data[ch]['good_sp'] != []:
					imind = np.array(data[ch]['g_sp_im_ind'])
					t_rel = np.array(data[ch]['g_sp_t_rel'])
					sp_ind = np.flatnonzero(imind == c_im_ind+first_img_ind)
					if sp_ind != [] and im == []:
						imname = data[ch]['g_sp_im_id'][sp_ind[0]]
						print imname
						im = pl.imread(imdir+imname+'.png')
					ax.plot(t_rel[sp_ind]/1000,ch_num*np.ones(np.shape(sp_ind)),'k.')
					ch_num += 1
			aim = pl.subplot(2,5,c_im_ind)
			aim.set_axis_off()
			aim.imshow(im,extent=[-25,175,0,200])
			aim.set_xlim([-50,200])
			aim.set_ylim([-50,200])
			ax.set_xlim([-50,200])
			ax.set_ylim([-1,43])
			#im_a = pl.axes([0.2,0.8,0.25,0.25])
			#ax.set_xticks([-100,0,100])
			if c_im_ind != 1:
				adjust_spines(ax,['bottom'])
			else:
				adjust_spines(ax,['left','bottom'])
				ax.set_xlabel('Time (ms)')
				ax.set_ylabel('Channel #')
			ax.set_xticks([0,100])
			ax.set_xticklabels(['0','100'])
			ax.grid(color = 'r',linestyle='-')
			ax.xaxis.grid(True)
			ax.yaxis.grid(False)
		pl.savefig('rasters_imind_'+str(first_img_ind)+'.pdf')
		
		
def plot_raster_byimg(data,imdir,channels):
	for im_num in [3,41,42,46,56,75,81]:#range(1,106,1): 
		im_id = 'Nat300_'+str(im_num)
		temp_id = np.flatnonzero((np.array(data[17]['b_sp_im_id']) ==im_id))
		temp_ind = np.array(data[17]['b_sp_im_ind'])
		im_ind = np.unique(temp_ind[temp_id])
		pl.figure()
		pl.clf()
		pl_num = 1
		fig_num = 1 
		for p_ind,c_im_ind in enumerate(im_ind):
			ch_num = 1
			ax = pl.subplot(2,5,pl_num+5)
			im = []
			for ch in channels:#data.keys():
				if data[ch]['good_sp'] != []:
					imind = np.array(data[ch]['g_sp_im_ind'])
					t_rel = np.array(data[ch]['g_sp_t_rel'])
					sp_ind = np.flatnonzero(imind == c_im_ind)
					if sp_ind != [] and im == []:
						imname = data[ch]['g_sp_im_id'][sp_ind[0]]
						im = pl.imread(imdir+imname+'.png')
					ax.plot(t_rel[sp_ind]/1000,ch_num*np.ones(np.shape(sp_ind)),'k.')
					ch_num += 1
			aim = pl.subplot(2,5,pl_num)
			aim.set_axis_off()
			aim.imshow(im,extent=[-25,175,0,200])
			aim.set_xlim([-50,200])
			aim.set_ylim([-50,200])
			ax.set_xlim([-50,200])
			ax.set_ylim([-1,43])
			#im_a = pl.axes([0.2,0.8,0.25,0.25])
			#ax.set_xticks([-100,0,100])
			if pl_num != 1:
				adjust_spines(ax,['bottom'])
			else:
				adjust_spines(ax,['left','bottom'])
				ax.set_xlabel('Time (ms)')
				ax.set_ylabel('Channel #')
			ax.set_xticks([0,100])
			ax.set_xticklabels(['0','100'])
			ax.grid(color = 'r',linestyle='-')
			ax.xaxis.grid(True)
			ax.yaxis.grid(False)
			pl_num += 1
			if pl_num > 5:
				pl.savefig('rasters_imd_'+str(im_id)+'_'+str(fig_num)+'.pdf')
				pl.clf()
				pl_num = 1
				fig_num += 1
		pl.close()
		

def plot_fp(data):
	im_id = np.array(data[17]['g_sp_im_id'])
	t_rel = np.array(data[17]['g_sp_t_rel'])
	fp = []
	for im_num in range(1,106,1):
		im_name = 'Nat300_'+str(im_num)
		sp_ind = np.flatnonzero(im_id==im_name)
		fp_a = len(np.flatnonzero((t_rel[sp_ind] > 75000 and t_rel[sp_ind] < 175000)))
		fp.append(len(sp_ind))
	pl.plot(fp)
	pl.plot(fp_a,'r')
	pl.show()
def main():
	data = pk.load(open('temp_data.pk','r'))
	imdir = '/Users/labuser/Desktop/ChaboSetup/NATURAL300/'
	#plot_crude_fp(data)
	#ch_clustered = cluster_neurons(data)
	#plot_crude_fp(data,ch_clustered)
	#plot_raster_byimg(data,imdir,ch_clustered)
	#plot_fp(data)
	get_fp(data)

			
	
	print 'done'
	pass


if __name__ == '__main__':
	main()

