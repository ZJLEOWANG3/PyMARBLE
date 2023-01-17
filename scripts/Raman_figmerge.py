#!/usr/bin/env python3
#############################################
# This script is to save and combine figures for cluster section with 
# different discussions
#############################################

from load import *

def savefig(figs, parentp, namelist):
    """
    to save figure
    input 
    figs : a list of python figure object
    parentp : parent path to save figure
    namelist : name list for the figures
    """
    for i in range(len(namelist)):
        figs[i].savefig(parentp+namelist[i],bbox_inches='tight',dpi=600)

def savedata(data, parentp, namelist):
    """
    to save data with correspondent figure
    input 
    data : a list of dataframe
    parentp : parent path to save df
    namelist : name list for the csv
    """
    for i in range(len(data)):
        try:
            data[i].to_csv(parentp+namelist[i])
        except AttributeError:
            np.save((parentp+namelist[i]).replace("csv","npy"),data[i])

def combinefig(im_list,outpp,outp,VOH,rowcol = 1, title = None,size=(800,200),ratio=None):
    """
    to combine figure
    outpp : parent output path
    outp : typically the name of the output file. e.g., a.png
    VOH : vertical "V" or horizontal "H"
    title : figure title list [ title position, title text, color, text size], if none, then add nothing
    rowcol : how many rows for vertical merge and how many cols for horizontal merge
    size : single size of each figure
    ratio: if None, each figure is renormalized into the same size
            but if a list, each figure is rescaled by the scalar pairs (width, height) in the list. e.g., (1,1) denotes unchange, (1,2) to double the height

    e.g., 
    im_list = [ Image.open(parentp+name) for name in namelist ] # image list
    Raman_figmerge.combinefig(im_list,outpp,outp,"H")
    """
    #transform images to same size
    ims = []
    for i in range(len(im_list)):
        if ratio ==None:
            temp_img = im_list[i].resize(size,Image.BILINEAR)
        else:
            temp_img = im_list[i].resize( (int(size[0]*ratio[i][0]),int(size[1]*ratio[i][1]) ),Image.BILINEAR)
        ims.append(temp_img)
        
    #create new blank figure
    t_width, t_height = 0, 0 # total width and height
    width,height = ims[0].size
    if VOH == "V":
        if ratio==None:
            t_width, t_height = width*rowcol,height*math.ceil(len(ims)/rowcol)
        else:
            for (w,h) in ratio:
                t_height += height * h
            t_width = max(ratio, key=operator.itemgetter(0))[0] * width * rowcol
    elif VOH == "H":
        if ratio==None:
            t_width, t_height = width*math.ceil(len(ims)/rowcol),height * rowcol
        else:
            for (w,h) in ratio:
                t_width += width * w
            t_height = max(ratio, key=operator.itemgetter(1))[1] * height * rowcol
    result = Image.new(ims[0].mode,(int(t_width),int(t_height)))


    #add each figure to the new one
    threshold = math.ceil(len(ims)/rowcol) # which rowcol to plot
    for i, im in enumerate(ims):
        if VOH == "V":
            if ratio==None:
                if i< threshold:
                    result.paste(im,box=(0,i*height))
                else:
                    result.paste(im,box=(int(t_width/rowcol),(i-threshold)*height))
            else:
                # start position of vertical = height*sum(h for _,h in ratio[:i])
                # start position of horizontal = -(w-t_width)/2 for w,_ in ratio[i]
                if i < threshold:
                    result.paste(im,box=( int((t_width/rowcol - ratio[i][0])/2),
                                        int(height*sum(h for _,h in ratio[:i])) ))
                else:
                    result.paste(im,box=( int((t_width/rowcol + ratio[i][0])/2),
                                        int(height*sum(h for _,h in ratio[ threshold+1:i])) ))
        elif VOH == "H":
            if ratio == None:
                if i< threshold:
                    result.paste(im,box=(i*width, 0))
                else:
                    result.paste(im,box=((i-threshold)*width, int(t_height/rowcol)))
            else:
                if i< threshold:
                    result.paste(im,box=( int(width*sum(w for w,_ in ratio[:i])),
                                        int((t_height/rowcol-ratio[i][1]*height)/2)  ))
                else:
                    result.paste(im,box=( int(width*sum(w for w,_ in ratio[threshold+1:i])),
                                        int((t_height/rowcol + ratio[i][1]*height)/2)  ))

    # save 
    result.save(os.path.join(outpp,outp))
    # add title
    """
    ImageDraw.Draw(
    result  # Image
    ).text(
    title[0],  # Coordinates
    title[1],  # Text
    title[2],  # Color
    fontsize = title[3]
    )
    result.save(outpp+outp)
    """



