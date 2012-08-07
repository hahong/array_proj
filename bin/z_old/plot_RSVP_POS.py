#!/usr/bin/env python
# encoding: utf-8
"""
plot_RSVP_POS.py

Created by najiblocal on 2010-09-28.
Copyright (c) 2010 __MyCompanyName__. All rights reserved.
"""

import sys
import getopt
import cPickle as pk
import pylab as pl
import numpy as np
from scipy.stats.stats import pearsonr

help_message = '''
The help message goes here.
'''


class Usage(Exception):
	def __init__(self, msg):
		self.msg = msg


def adjust_spines(ax,spines):
    for loc, spine in ax.spines.iteritems():
        if loc in spines:
            spine.set_position(('outward',2)) # outward by 10 points
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

def data_to_ch(data):
	ch = {}
	for ch_ind in range(1,97):
		ch[ch_ind] = {}
		ch[ch_ind]['bl'] = data [ch_ind]['blanks']
		ch[ch_ind]['bl_mu'] = pl.mean(ch[ch_ind]['bl'])
		ch[ch_ind]['bl_sem'] =  pl.std(ch[ch_ind]['bl'])/pl.sqrt(len(ch[ch_ind]['bl']))
		for ind in sorted(data[ch_ind].keys()):
			if ind != 'blanks':
				k = ind[0]
				if k not in ch[ch_ind]:
					ch[ch_ind][k] = {}
					ch[ch_ind][k]['fr'] = []
					ch[ch_ind][k]['fr_mu'] = []
					ch[ch_ind][k]['fr_sem'] = []
					ch[ch_ind][k]['pos_y'] = []
					ch[ch_ind][k]['dprime']=[]
				ch[ch_ind][k]['fr'].append(data[ch_ind][ind]['on'])
				ch[ch_ind][k]['fr_mu'].append(pl.mean(data[ch_ind][ind]['on']))
				ch[ch_ind][k]['fr_sem'].append(pl.std(data[ch_ind][ind]['on'])/pl.sqrt(len(data[1][ind]['on'])))
				ch[ch_ind][k]['pos_y'].append(ind[2])
				#print ch[ch_ind][k]['pos_y']
				#print pl.std(data[ch_ind][ind]['on'])
				ch[ch_ind][k]['dprime'].append((pl.mean(data[ch_ind][ind]['on'])-ch[ch_ind]['bl_mu'])/((pl.std(ch[ch_ind]['bl'])+pl.std(data[ch_ind][ind]['on']))/2))
		#print ch[ch_ind]['OSImage_5']['pos_y']
	return ch
	
def plot_all_channels_assignPN(ch1,ch2,sig_ch_ind,images,swap_pos,fn1,fn2):
	im_a = images[0]
	im_b = images[1]
	fig_pre = pl.figure(1)
	fig_post = pl.figure(2)
	fig_hist = pl.figure(3)
	btstrpitr = 1
	deltapminusn = np.zeros([btstrpitr,3])
	print len(sig_ch_ind)
	if len(sig_ch_ind) > 24:
		num_sps = 24
	else: 
		num_sps = len(sig_ch_ind)
	swap_dPN = []
	n_swap_dPN = []
	p_ch = []
	fov_dPN = []
	for t_ind in range(num_sps):
		ch_ind = sig_ch_ind[t_ind]
		ax_pre = fig_pre.add_subplot(5,5,t_ind+1)
		ax_post = fig_post.add_subplot(5,5,t_ind+1)
		ax_hist = fig_hist.add_subplot(5,5,t_ind+1)
		for a in range(btstrpitr):
			im_a_1 = np.zeros(3)
			im_a_2 = np.zeros(3)
			im_b_1 = np.zeros(3)
			im_b_2 = np.zeros(3)
			for pos_ind in range(3):
				rand_ind_a = np.random.permutation(len(ch1[ch_ind][im_a]['fr'][pos_ind]))
				im_a_1[pos_ind] = np.mean(ch1[ch_ind][im_a]['fr'][pos_ind][rand_ind_a[:20]])
				im_a_2[pos_ind] = np.mean(ch1[ch_ind][im_a]['fr'][pos_ind][:])
				rand_ind_b = np.random.permutation(len(ch1[ch_ind][im_b]['fr'][pos_ind]))
				im_b_1[pos_ind] = np.mean(ch1[ch_ind][im_b]['fr'][pos_ind][rand_ind_b[:20]])
				im_b_2[pos_ind] = np.mean(ch1[ch_ind][im_b]['fr'][pos_ind][:])
				
			if im_a_1[1] > im_b_1[1]:
				im_p_fr = im_a_2
				im_n_fr = im_b_2
				im_p = im_a
				im_n = im_b				
			else:
				im_p_fr = im_b_2
				im_n_fr = im_a_2
				im_p = im_b
				im_n = im_a
			if im_p_fr[1]-im_n_fr[1] > 40:
				#rand_ind_blank = np.random.permutation(len(ch1[ch_ind]['bl']))
				#bl_fr = np.mean(ch1[ch_ind]['bl'][rand_ind_blank[:30]])
				#ax_pre.plot(ch1[ch_ind][im_a]['pos_y'],bl_fr*np.ones((3,1)).flatten(),'k')
				ax_pre.errorbar(ch1[ch_ind][im_a]['pos_y'],ch1[ch_ind][im_p]['fr_mu'],ch1[ch_ind][im_p]['fr_sem'],color = 'c',capsize= 0)
				ax_pre.errorbar(ch1[ch_ind][im_b]['pos_y'],ch1[ch_ind][im_n]['fr_mu'],ch1[ch_ind][im_n]['fr_sem'],color = 'r',capsize= 0)
				ax_post.errorbar(ch2[ch_ind][im_p]['pos_y'],ch2[ch_ind][im_p]['fr_mu'],ch2[ch_ind][im_p]['fr_sem'],color = 'c',capsize= 0)
				ax_post.errorbar(ch2[ch_ind][im_n]['pos_y'],ch2[ch_ind][im_n]['fr_mu'],ch2[ch_ind][im_n]['fr_sem'],color = 'r',capsize= 0)
				#ax_post.errorbar(ch2[ch_ind][im_n]['pos_y'],ch2[ch_ind]['bl_mu']*np.ones((3,1)).flatten(),yerr = ch2[ch_ind]['bl_sem'],capsize = 0,color = 'k')
				deltapminusn[a,:] = (np.array(ch2[ch_ind][im_p]['fr_mu'])-np.array(ch2[ch_ind][im_n]['fr_mu']))-(im_p_fr-im_n_fr)
				p_ch.append(ch_ind)
				if swap_pos == -3:
					ax_hist.hist(deltapminusn[:,0],bins = range(-100,100,5),facecolor='r',edgecolor = 'none',alpha = 0.25)
					ax_hist.hist(deltapminusn[:,2],bins = range(-100,100,5),facecolor='none',edgecolor='g')
					ax_hist.set_xlim([-50,50])
					swap_dPN.append(np.mean(deltapminusn[:,0]))
					n_swap_dPN.append(np.mean(deltapminusn[:,2]))
					fov_dPN.append(np.mean(deltapminusn[:,1]))	
				elif swap_pos == 3:
					ax_hist.hist(deltapminusn[:,2],bins = range(-100,100,5),facecolor='r',edgecolor = 'none',alpha = 0.25)
					ax_hist.hist(deltapminusn[:,0],bins = range(-100,100,5),facecolor='none',edgecolor='g')
					ax_hist.set_xlim([-50,50])
					swap_dPN.append(np.mean(deltapminusn[:,2]))
					n_swap_dPN.append(np.mean(deltapminusn[:,0]))
					fov_dPN.append(np.mean(deltapminusn[:,1]))	
			else:
				ax_pre.errorbar(ch1[ch_ind][im_a]['pos_y'],ch1[ch_ind][im_p]['fr_mu'],ch1[ch_ind][im_p]['fr_sem'],color = 'c',capsize= 0,alpha= 0.1)
				ax_pre.errorbar(ch1[ch_ind][im_b]['pos_y'],ch1[ch_ind][im_n]['fr_mu'],ch1[ch_ind][im_n]['fr_sem'],color = 'r',capsize= 0,alpha= 0.1)
				ax_post.errorbar(ch2[ch_ind][im_p]['pos_y'],ch2[ch_ind][im_p]['fr_mu'],ch2[ch_ind][im_p]['fr_sem'],color = 'c',capsize= 0,alpha= 0.1)
				ax_post.errorbar(ch2[ch_ind][im_n]['pos_y'],ch2[ch_ind][im_n]['fr_mu'],ch2[ch_ind][im_n]['fr_sem'],color = 'r',capsize= 0,alpha= 0.1)
		#fix graph
		#pl.errorbar(ch[ch_ind][im_a]['pos_y'],ch[ch_ind]['bl_mu']*np.ones((3,1)).flatten(),yerr = ch[ch_ind]['bl_sem'],capsize = 0,color = 'k')
		all_fr = ch1[ch_ind][im_a]['fr_mu']+ch1[ch_ind][im_b]['fr_mu']
		all_fr.append(ch1[ch_ind]['bl_mu'])
		all_sem = ch1[ch_ind][im_a]['fr_sem']+ch1[ch_ind][im_b]['fr_sem']
		all_sem.append(ch1[ch_ind]['bl_sem'])
		ax_pre.text(4,max(all_fr),str(sig_ch_ind[t_ind]))
		ax_pre.set_xlim([-5,5])
		ax_pre.set_xticks([-3,0,3])
		ax_pre.set_ylim([min(all_fr)-max(all_sem),max(all_fr)+max(all_sem)])
		all_fr_post = ch2[ch_ind][im_a]['fr_mu']+ch2[ch_ind][im_b]['fr_mu']
		all_fr_post.append(ch2[ch_ind]['bl_mu'])
		all_sem_post = ch2[ch_ind][im_a]['fr_sem']+ch2[ch_ind][im_b]['fr_sem']
		all_sem_post.append(ch2[ch_ind]['bl_sem'])
		ax_post.text(4,max(all_fr_post),str(sig_ch_ind[t_ind]))
		ax_post.plot([swap_pos,swap_pos],[min(all_fr_post)-max(all_sem_post),max(all_fr_post)+max(all_sem_post)],lw=4,color='r',alpha=0.25)
		ax_post.set_xlim([-5,5])
		ax_post.set_xticks([-3,0,3])
		ax_post.set_ylim([min(all_fr_post)-max(all_sem_post),max(all_fr_post)+max(all_sem_post)])
		adjust_spines(ax_pre,['left','bottom'])
		adjust_spines(ax_post,['left','bottom'])
	ax_sum = fig_hist.add_subplot(5,5,25)
	ax_sum.hist(swap_dPN,bins = range(-100,100,5),facecolor='r',edgecolor = 'none',alpha = 0.25)
	ax_sum.hist(n_swap_dPN,bins = range(-100,100,5),facecolor='none',edgecolor='g')
	ax_sum.set_xlim([-50,50])
	t = fn1.split('.')[0] + 'v' + fn2.split('.')[0].split('_')[-1]
	fig_pre.suptitle(t+'_pre')
	fig_post.suptitle(t+'_post')
	fig_hist.suptitle(t)
	fig_pre.savefig('./plots/'+t+'_pre.pdf')
	fig_post.savefig('./plots/'+t+'_post.pdf')
	fig_hist.savefig('./plots/'+t+'_dPN.pdf')
	pf = open('./out_njm/'+ t + '_out.pk','wb')
	pk.dump({'fov_dPN':fov_dPN,'s_dPN':swap_dPN,'ns_dPN':n_swap_dPN,'ch':p_ch},pf)
	pf.close()

def plot_average_channels(ch):
	avfr = {}
	avdp = {}
	for ind,ch_ind in enumerate(ch):
		#print ind,ch_ind
		for k in sorted(ch[ch_ind].keys()):
			if k not in avfr and 'bl' not in k:
				avfr[k] = pl.empty([96,3])
				avdp[k] = pl.empty([96,3])
			if 'bl' not in k:
				avfr[k][ind,:] = ch[ch_ind][k]['fr_mu']-ch[ch_ind]['bl_mu']
				avdp[k][ind,:] = ch[ch_ind][k]['dprime']
	pl.subplot(3,1,1)
	pl.hist(avfr['OSImage_5'][:,0]-avfr['OSImage_45'][:,0],bins=range(-30,30,5))
	pl.subplot(3,1,2)
	pl.hist(avfr['OSImage_5'][:,1]-avfr['OSImage_45'][:,1],bins=range(-30,30,5))
	pl.subplot(3,1,3)
	pl.hist(avfr['OSImage_5'][:,2]-avfr['OSImage_45'][:,2],bins=range(-30,30,5))

def get_pearson_corr(ch1,ch2,images,images_l,fn1,fn2):
	#computer correlation for non mainpulated images.
	ch1fp = {}
	ch2fp = {}
	ch1l = {}
	ch2l = {}
	for ch_ind in range(1,97):
		ch1fp[ch_ind] = []
		ch2fp[ch_ind] = []
		ch1l[ch_ind] = []
		ch2l[ch_ind] = []
		#ch1fp[ch_ind].append(ch1[ch_ind]['bl_mu'])
		#ch2fp[ch_ind].append(ch2[ch_ind]['bl_mu'])
		for imname in images:
			ch1fp[ch_ind].append(ch1[ch_ind][imname]['fr_mu'])#-ch1[ch_ind]['bl_mu'])
			ch2fp[ch_ind].append(ch2[ch_ind][imname]['fr_mu'])#-ch1[ch_ind]['bl_mu'])
		for imname in images_l:
			ch1l[ch_ind].append(ch1[ch_ind][imname]['fr_mu'])
			ch2l[ch_ind].append(ch2[ch_ind][imname]['fr_mu'])
	chrho = []
	for ch_ind in range(1,97):
		#print ch1fp[ch_ind]
		#print np.array(ch2fp[ch_ind]).flatten()
		chrho.append(pearsonr(np.array(ch1fp[ch_ind]).flatten(),np.array(ch2fp[ch_ind]).flatten())[0])
	#print np.array(chrho)
	#print pl.shape(np.array(chrho))
	fig_corr_a = pl.figure(5)
	s_p_a = fig_corr_a.add_subplot(1,2,1)
	s_p_a.imshow(np.reshape(np.array(chrho),(8,12)),norm = None,vmin = 0, vmax = 1,cmap = pl.cm.gray,interpolation = None)	
	s_p_b = fig_corr_a.add_subplot(1,2,2)
	s_p_b.hist(np.array(chrho),bins = np.arange(-1.1,1.1,0.1))
	sig_channels = pl.nonzero(np.array(chrho) > 0.85)[0]+1
	fig_corr_b = pl.figure(6)
	for n_plot in range(len(sig_channels)):
		s_p_temp = fig_corr_b.add_subplot(5,5,n_plot+1)
		s_p_temp.plot(np.array(ch1fp[sig_channels[n_plot]]).flatten(),np.array(ch2fp[sig_channels[n_plot]]).flatten(),'k.')
		s_p_temp.plot(np.array(ch1l[sig_channels[n_plot]]).flatten(),np.array(ch2l[sig_channels[n_plot]]).flatten(),'r.')
		s_p_temp.plot(ch1[sig_channels[n_plot]]['bl_mu'],ch2[sig_channels[n_plot]]['bl_mu'],'kx')
		s_p_temp.plot([0,100],[0,100],'k-')
	t = fn1.split('.')[0] + 'v' + fn2.split('.')[0].split('_')[-1]
	fig_corr_a.savefig('./plots/'+t+'corr_a.pdf')
	fig_corr_b.savefig('./plots/'+t+'corr_b.pdf')
	return sig_channels
	

def main(argv=sys.argv):
	print len(argv)
	if len(argv) !=3 and len(argv)!=1:
		print argv[0] + " <input1> <input2>"
		return 2
	fn1 = []
	fn2 = []
	if len(argv) == 3:
		fn1.append(argv[1])
		fn2.append(argv[2])
	
	if len(argv) == 1:
		for line in open('pre.txt'):
			if '\n'in line:fn1.append(line[:-1])
			else:fn1.append(line)
		for line in open('post.txt'):
			if '\n' in line:fn2.append(line[:-1])
			else:fn2.append(line)
	print fn1
	print fn2			
	for f_index in range(len(fn1)):
		filename1 = fn1[f_index].split('/')[-1]
		filename2 = fn2[f_index].split('/')[-1]
		data1 = pk.load(open('./data_postproc/'+filename1))
		ch1 = data_to_ch(data1)
		data2 = pk.load(open('./data_postproc/'+filename2))
		ch2 = data_to_ch(data2)
		print filename1.split('_')[1]
		if filename1.split('_')[1] <= '20101001':
			images_fp = ['OSImage_16','OSImage_32','OSImage_51']
			images_l = ['OSImage_45','OSImage_5']
			swap_pos = -3
		elif filename1.split('_')[1] > '20101001':
			images_fp = ['OSImage_45','OSImage_5','OSImage_51']
			images_l = ['OSImage_16','OSImage_32']
			swap_pos = 3
		sig_ch_ind = get_pearson_corr(ch1,ch2,images_fp,images_l,filename1,filename2)
		print sig_ch_ind
		plot_all_channels_assignPN(ch1,ch2,sig_ch_ind,images_l,swap_pos,filename1,filename2)
		#pl.show()
		pl.close('all')	


if __name__ == "__main__":
	sys.exit(main())
