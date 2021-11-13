

def update_wirc(img):
    
    header = img[0].header
    
    header["FILTER"] = header["AFT"].split("__")[0]
    
    img[0].header = header
    
    return img