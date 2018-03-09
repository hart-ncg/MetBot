# DICTIONARY OF FILTER VALUES
blobfilters={}
blobfilters['SAcloudband']={
'noaa-olr-0-0': {'area': (4,3000), 'angle': (-90., -5.0), 'latextent': (-20,-40), 'ROI': (7.5, 100,-15,-40),'angleROI': (7.5, 100,-23,-33), 'stretch': (100,320), 'thresh': (230,'low'), 'text': "Cloudband" },
'noaa-olr-0-all': {'area': (4,3000), 'angle': False, 'latextent': False, 'ROI': (7.5, 100,-15,-40), 'stretch': (100,320), 'thresh': (230,'low'), 'text': "Cloudband" },
'noaa-olr-0-full': {'area': (4,3000), 'angle': False, 'latextent': False, 'ROI': False, 'stretch': (100,320), 'thresh': (230,'low'), 'text': "Cloudband" },

'ncdr-olr-0-0': {'area': (4,3000), 'angle': (-90., -5.0), 'latextent': (-20,-40), 'ROI': (7.5, 100,-15,-40),'angleROI': (7.5, 100,-23,-33), 'stretch': (100,320), 'thresh': (240,'low'), 'text': "Cloudband" },
'ncdr-olr-0-all': {'area': (4,3000), 'angle': False, 'latextent': False, 'ROI': (7.5, 100,-15,-40), 'stretch': (100,320), 'thresh': (240,'low'), 'text': "Cloudband" },
'ncdr-olr-0-full': {'area': (4,3000), 'angle': False, 'latextent': False, 'ROI': False, 'stretch': (100,320), 'thresh': (240,'low'), 'text': "Cloudband" },

'ncep2-olr-0-0': {'area': (4,3000), 'angle': (-90., -5.0), 'latextent': (-20,-40), 'ROI': (7.5, 100,-15,-40),'angleROI': (7.5, 100,-23,-33), 'stretch': (100,320), 'thresh': (245,'low'), 'text': "Cloudband" },
'ncep2-olr-0-all': {'area': (4,3000), 'angle': False, 'latextent': False, 'ROI': (7.5, 100,-15,-40), 'stretch': (100,320), 'thresh': (245,'low'), 'text': "Cloudband" },
'ncep2-olr-0-full': {'area': (4,3000), 'angle': False, 'latextent': False, 'ROI': False, 'stretch': (100,320), 'thresh': (245,'low'), 'text': "Cloudband" },

'erain-olr-0-0': {'area': (4,3000), 'angle': (-90., -5.0), 'latextent': (-20,-40), 'ROI': (7.5, 100,-15,-40),'angleROI': (7.5, 100,-23,-33), 'stretch': (100,320), 'thresh': (245,'low'), 'text': "Cloudband" },
'erain-olr-0-all': {'area': (4,3000), 'angle': False, 'latextent': False, 'ROI': (7.5, 100,-15,-40), 'stretch': (100,320), 'thresh': (245,'low'), 'text': "Cloudband" },
'erain-olr-0-full': {'area': (4,3000), 'angle': False, 'latextent': False, 'ROI': False, 'stretch': (100,320), 'thresh': (245,'low'), 'text': "Cloudband" },

'ncep2-hgt-850-del2': {'area': (1000, 1000000), 'angle': False, 'latextent': False, 'ROI': False, 'stretch': (-3,4), 'thresh': (1.0,'high'), 'text': "low"},
'ncep2-hgt-500-del2': {'area': (1000, 1000000), 'angle': False, 'latextent': False, 'ROI': False, 'stretch': (-3.5,4.5), 'thresh': (1,'high'), 'text': "depression"},
'ncep2-hgt-250-del2': {'area': (1000, 1000000), 'angle': False, 'latextent': False, 'ROI': False, 'stretch': (-6,7), 'thresh': (1,'high'), 'text': "trough"},
'ncep2-uvmag-500-del2': {'area': (1000, 1000000), 'angle': False, 'latextent': False, 'ROI': False, 'stretch': (-2.5,2.5), 'thresh': (-1,'low'), 'text': "jet"},
'ncep2-uvmag-250-del2': {'area': (1000, 1000000), 'angle': False, 'latextent': False, 'ROI': False, 'stretch': (-2.5,2.5), 'thresh': (-1,'low'), 'text': "jet"},
'ncep2-strmfp-200-0': {'area': (15, 5000), 'angle': False, 'latextent': False, 'ROI': False, 'stretch': (-2.5,2.5), 'thresh': (0.2e7,'high'), 'text': "eddy"},
'ncep2-strmfn-200-0': {'area': (15, 5000), 'angle': False, 'latextent': False, 'ROI': False, 'stretch': (-2.5,2.5), 'thresh': (-0.2e7,'low'), 'text': "eddy"},
'ncep2-uvmag-850-del2': {'area': (1000, 1000000), 'angle': False, 'latextent': False, 'ROI': False, 'stretch': (-1.5,1.5), 'thresh': (-0.5,'low'), 'text': "jet"},
'ncep2-quvmag-500-del2': {'area': (1000, 1000000), 'angle': False, 'latextent': False, 'ROI': False, 'stretch': (-2.5,2.5), 'thresh': (-1,'low'), 'text': "moistjet"},
'ncep2-quvmag-700-del2': {'area': (1000, 1000000), 'angle': False, 'latextent': False, 'ROI': False, 'stretch': (-2.5,2.5), 'thresh': (-1,'low'), 'text': "moistjet"},
'ncep2-quvmag-850-del2': {'area': (1000, 1000000), 'angle': False, 'latextent': False, 'ROI': False, 'stretch': (-1.5,1.5), 'thresh': (-0.5,'low'), 'text': "moistjet"},

'ncep-olr-0-0': {'area': (4,3000), 'angle': (-90., -5.0), 'latextent': (-20,-40), 'ROI': (7.5, 100,-15,-40),'angleROI': (7.5, 100,-23,-33), 'stretch': (100,320), 'thresh': (245,'low'), 'text': "Cloudband" },
'ncep-olr-0-all': {'area': (4,3000), 'angle': False, 'latextent': False, 'ROI': (7.5, 100,-15,-40), 'stretch': (100,320), 'thresh': (245,'low'), 'text': "Cloudband" },
'ncep-olr-0-full': {'area': (4,3000), 'angle': False, 'latextent': False, 'ROI': False, 'stretch': (100,320), 'thresh': (245,'low'), 'text': "Cloudband" },

'era-olr-0-0': {'area': (4,3000), 'angle': (-90., -5.0), 'latextent': (-20,-40), 'ROI': (7.5, 100,-15,-40),'angleROI': (7.5, 100,-23,-33), 'stretch': (100,320), 'thresh': (245,'low'), 'text': "Cloudband" },
'era-olr-0-all': {'area': (4,3000), 'angle': False, 'latextent': False, 'ROI': (7.5, 100,-15,-40), 'stretch': (100,320), 'thresh': (245,'low'), 'text': "Cloudband" },
'era-olr-0-full': {'area': (4,3000), 'angle': False, 'latextent': False, 'ROI': False, 'stretch': (100,320), 'thresh': (245,'low'), 'text': "Cloudband" },

'20cr-olr-0-0': {'area': (4,3000), 'angle': (-90., -5.0), 'latextent': (-20,-40), 'ROI': (7.5, 100,-15,-40),'angleROI': (7.5, 100,-23,-33), 'stretch': (100,320), 'thresh': (240,'low'), 'text': "Cloudband" },
'20cr-olr-0-all': {'area': (4,3000), 'angle': False, 'latextent': False, 'ROI': (7.5, 100,-15,-40), 'stretch': (100,320), 'thresh': (240,'low'), 'text': "Cloudband" },
'20cr-olr-0-full': {'area': (4,3000), 'angle': False, 'latextent': False, 'ROI': False, 'stretch': (100,320), 'thresh': (240,'low'), 'text': "Cloudband" },

'um-olr-0-0': {'area': (4,3000), 'angle': (-90., -5.0), 'latextent': (-20,-40), 'ROI': (7.5, 100,-15,-40),'angleROI': (7.5, 100,-23,-33), 'stretch': (100,320), 'thresh': (245,'low'), 'text': "Cloudband" },
'um-olr-0-all': {'area': (4,3000), 'angle': False, 'latextent': False, 'ROI': (7.5, 100,-15,-40), 'stretch': (100,320), 'thresh': (245,'low'), 'text': "Cloudband" },
'um-olr-0-full': {'area': (4,3000), 'angle': False, 'latextent': False, 'ROI': False, 'stretch': (100,320), 'thresh': (245,'low'), 'text': "Cloudband" },

'cmip5-olr-0-0': {'area': (4,3000), 'angle': (-90., -5.0), 'latextent': (-20,-40), 'ROI': (7.5, 100,-15,-40),'angleROI': (7.5, 100,-23,-33), 'stretch': (100,320), 'thresh': (245,'low'), 'text': "Cloudband" },
'cmip5-olr-0-all': {'area': (4,3000), 'angle': False, 'latextent': False, 'ROI': (7.5, 100,-15,-40), 'stretch': (100,320), 'thresh': (245,'low'), 'text': "Cloud" },
'cmip5-olr-0-full': {'area': (4,3000), 'angle': False, 'latextent': False, 'ROI': False, 'stretch': (100,320), 'thresh': (245,'low'), 'text': "Cloud" },

'cfsr-olr-0-0': {'area': (5000, 1000000), 'angle': (5.0, 90.0), 'latextent': (-23,-35), 'ROI': (7.5, 80,-23,-40),'angleROI': (7.5, 80,-23,-33), 'stretch': (100,320), 'thresh': (245,'low'), 'text': "Cloudband" },
'cfsr-olr-0-all': {'area': (5000, 1000000), 'angle': False, 'latextent': False, 'ROI': (7.5, 80,-23,-40), 'stretch': (100,320), 'thresh': (245,'low'), 'text': "Cloudband" },
'cfsr-olr-0-full': {'area': (5000, 1000000), 'angle': False, 'latextent': False, 'ROI': False, 'stretch': (100,320), 'thresh': (245,'low'), 'text': "Cloudband" },
'cfsr-olr-0-d': {'area': (5000, 1000000), 'angle': (5.0, 90.0), 'latextent': (-23,-35), 'ROI': (7.5, 80,-23,-40), 'stretch': (100,320), 'thresh': (245,'low'), 'text': "Cloudband" },
'cfsr-olr-0-dall': {'area': (5000, 1000000), 'angle': False, 'latextent': False, 'ROI': (7.5, 80,-23,-40), 'stretch': (100,320), 'thresh': (245,'low'), 'text': "Cloudband" },
'cfsr-olr-0-dfull': {'area': (5000, 1000000), 'angle': False, 'latextent': False, 'ROI': False, 'stretch': (100,320), 'thresh': (245,'low'), 'text': "Cloudband" },
'cfsr-hgt-850-del2': {'area': (1000, 1000000), 'angle': False, 'latextent': False, 'ROI': False, 'stretch': (-3,4), 'thresh': (1.0,'high'), 'text': "low"},
'cfsr-hgt-500-del2': {'area': (1000, 1000000), 'angle': False, 'latextent': False, 'ROI': False, 'stretch': (-3.5,4.5), 'thresh': (1,'high'), 'text': "depression"},
'cfsr-hgt-250-del2': {'area': (1000, 1000000), 'angle': False, 'latextent': False, 'ROI': False, 'stretch': (-6,7), 'thresh': (1,'high'), 'text': "trough"},
'cfsr-uvmag-500-del2': {'area': (1000, 1000000), 'angle': False, 'latextent': False, 'ROI': False, 'stretch': (-2.5,2.5), 'thresh': (-1,'low'), 'text': "jet"},
'cfsr-uvmag-250-del2': {'area': (1000, 1000000), 'angle': False, 'latextent': False, 'ROI': False, 'stretch': (-2.5,2.5), 'thresh': (-1,'low'), 'text': "jet"},
'cfsr-uvmag-850-del2': {'area': (1000, 1000000), 'angle': False, 'latextent': False, 'ROI': False, 'stretch': (-1.5,1.5), 'thresh': (-0.5,'low'), 'text': "jet"},
'cfsr-quvmag-500-del2': {'area': (1000, 1000000), 'angle': False, 'latextent': False, 'ROI': False, 'stretch': (-2.5,2.5), 'thresh': (-1,'low'), 'text': "moistjet"},
'cfsr-quvmag-700-del2': {'area': (1000, 1000000), 'angle': False, 'latextent': False, 'ROI': False, 'stretch': (-2.5,2.5), 'thresh': (-1,'low'), 'text': "moistjet"},
'cfsr-quvmag-850-del2': {'area': (1000, 1000000), 'angle': False, 'latextent': False, 'ROI': False, 'stretch': (-1.5,1.5), 'thresh': (-0.5,'low'), 'text': "moistjet"},

'hadam3p-olr-0-0': {'area': (5000, 1000000), 'angle': (5.0, 90.0), 'latextent': (-20,-40), 'ROI': (7.5, 78.75,-15,-40),'angleROI': (7.5, 78.75,-22.5,-33.75), 'stretch': (100,320), 'thresh': (245,'low'), 'text': "Cloudband" },
'hadam3p-olr-0-all': {'area': (5000, 1000000), 'angle': False, 'latextent': False, 'ROI': (7.5, 78.75,-15,-40), 'stretch': (100,320), 'thresh': (245,'low'), 'text': "Cloudband" },
'hadam3p-olr-0-full': {'area': (5000, 1000000), 'angle': False, 'latextent': False, 'ROI': False, 'stretch': (100,320), 'thresh': (245,'low'), 'text': "Cloudband" },
'hadam3p-hgt-850-del2': {'area': (1000, 1000000), 'angle': False, 'latextent': False, 'ROI': False, 'stretch': (-3,4), 'thresh': (1.0,'high'), 'text': "low"},
'hadam3p-hgt-500-del2': {'area': (1000, 1000000), 'angle': False, 'latextent': False, 'ROI': False, 'stretch': (-3.5,4.5), 'thresh': (1,'high'), 'text': "depression"},
'hadam3p-hgt-250-del2': {'area': (1000, 1000000), 'angle': False, 'latextent': False, 'ROI': False, 'stretch': (-6,7), 'thresh': (1,'high'), 'text': "trough"},
'hadam3p-uvmag-250-del2': {'area': (1000, 1000000), 'angle': False, 'latextent': False, 'ROI': False, 'stretch': (-2.5,2.5), 'thresh': (-1,'low'), 'text': "jet"},
'hadam3p-uvmag-850-del2': {'area': (1000, 1000000), 'angle': False, 'latextent': False, 'ROI': False, 'stretch': (-1.5,1.5), 'thresh': (-0.5,'low'), 'text': "jet"},
'hadam3p-quvmag-700-del2': {'area': (1000, 1000000), 'angle': False, 'latextent': False, 'ROI': False, 'stretch': (-2.5,2.5), 'thresh': (-1,'low'), 'text': "moistjet"},
'hadam3p-quvmag-850-del2': {'area': (1000, 1000000), 'angle': False, 'latextent': False, 'ROI': False, 'stretch': (-1.5,1.5), 'thresh': (-0.5,'low'), 'text': "moistjet"},

'cp4-olr-0-0': {'area': (4,3000), 'angle': (-90., -5.0), 'latextent': (-21,-39), 'ROI': (7.5, 56,-20,-40),'angleROI': (7.5, 56,-23,-33), 'stretch': (100,320), 'thresh': (245,'low'), 'text': "Cloudband" },
'cp4-olr-0-all': {'area': (4,3000), 'angle': False, 'latextent': False, 'ROI': (7.5, 56,-20,-40), 'stretch': (100,320), 'thresh': (245,'low'), 'text': "Cloudband" },
'cp4-olr-0-full': {'area': (4,3000), 'angle': False, 'latextent': False, 'ROI': False, 'stretch': (100,320), 'thresh': (245,'low'), 'text': "Cloudband" },

}

blobfilters['SH60cloudband']={
'noaa-olr-0-0': {'area': (5000, 1000000), 'angle': (-90., -5.0),
                 'latextent': (-20,-40), 'ROI': (0, 360,-15,-40),
                 'angleROI': (0, 360,-23,-33), 'stretch': (100,320),
                 'thresh': (230,'low'), 'text': "Cloudband" },
'noaa-olr-0-all': {'area': (5000, 1000000), 'angle': False, 
                   'latextent': False, 'ROI': (0, 360,-15,-40),
                   'stretch': (100,320), 'thresh': (230,'low'), 
                   'text': "Cloudband" },
'noaa-olr-0-full': {'area': (5000, 1000000), 'angle': False, 
                    'latextent': False, 'ROI': False, 'stretch': (100,320),
                    'thresh': (230,'low'), 'text': "Cloudband" },
'ncdr-olr-0-0': {'area': (4,3000), 'angle': (-90., -5.0), 'latextent': (-20,-40), 'ROI': (0, 360,-15,-40),'angleROI': (0, 360,-23,-33), 'stretch': (100,320), 'thresh': (240,'low'), 'text': "Cloudband" },
'ncdr-olr-0-all': {'area': (4,3000), 'angle': False, 'latextent': False, 'ROI': (0, 360,-15,-40), 'stretch': (100,320), 'thresh': (240,'low'), 'text': "Cloudband" },
'ncdr-olr-0-full': {'area': (4,3000), 'angle': False, 'latextent': False, 'ROI': False, 'stretch': (100,320), 'thresh': (240,'low'), 'text': "Cloudband" }}

blobfilters['NAcloudband']={
'noaa-olr-0-0': {'area': (5000, 1000000), 'angle': (-90.0, -5.0), 'latextent': (35,23), 'ROI': (280.0,357.5,40,20), 'stretch': (100,320), 'thresh': (230,'low'), 'text': "Cloudband" },
'noaa-olr-0-all': {'area': (5000, 1000000), 'angle': False, 'latextent': False, 'ROI': (7.5, 80,-23,-40), 'stretch': (100,320), 'thresh': (230,'low'), 'text': "Cloudband" },
'noaa-olr-0-full': {'area': (5000, 1000000), 'angle': False, 'latextent': False, 'ROI': False, 'stretch': (100,320), 'thresh': (230,'low'), 'text': "Cloudband" },
'ncep2-hgt-850-del2': {'area': (1000, 1000000), 'angle': False, 'latextent': False, 'ROI': False, 'stretch': (-3,4), 'thresh': (1.0,'high'), 'text': "low"},
'ncep2-hgt-500-del2': {'area': (1000, 1000000), 'angle': False, 'latextent': False, 'ROI': False, 'stretch': (-3.5,4.5), 'thresh': (1,'high'), 'text': "depression"},
'ncep2-hgt-250-del2': {'area': (1000, 1000000), 'angle': False, 'latextent': False, 'ROI': False, 'stretch': (-6,7), 'thresh': (1,'high'), 'text': "trough"},
'ncep2-uvmag-250-del2': {'area': (1000, 1000000), 'angle': False, 'latextent': False, 'ROI': False, 'stretch': (-2.5,2.5), 'thresh': (-1,'low'), 'text': "jet"},
'ncep2-uvmag-850-del2': {'area': (1000, 1000000), 'angle': False, 'latextent': False, 'ROI': False, 'stretch': (-1.5,1.5), 'thresh': (-0.5,'low'), 'text': "jet"},
'ncep2-quvmag-700-del2': {'area': (1000, 1000000), 'angle': False, 'latextent': False, 'ROI': False, 'stretch': (-2.5,2.5), 'thresh': (-1,'low'), 'text': "moistjet"},
'ncep2-quvmag-850-del2': {'area': (1000, 1000000), 'angle': False, 'latextent': False, 'ROI': False, 'stretch': (-1.5,1.5), 'thresh': (-0.5,'low'), 'text': "moistjet"},
}

blobfilters['NPcloudband']={
'noaa-olr-0-0': {'area': (5000, 1000000), 'angle': (-90.0, -5.0), 'latextent': (35,23), 'ROI': (180.0,260.0,40,20), 'stretch': (100,320), 'thresh': (230,'low'), 'text': "Cloudband" },
'noaa-olr-0-all': {'area': (5000, 1000000), 'angle': False, 'latextent': False, 'ROI': (7.5, 80,-23,-40), 'stretch': (100,320), 'thresh': (230,'low'), 'text': "Cloudband" },
'noaa-olr-0-full': {'area': (5000, 1000000), 'angle': False, 'latextent': False, 'ROI': False, 'stretch': (100,320), 'thresh': (230,'low'), 'text': "Cloudband" },
'ncep2-hgt-850-del2': {'area': (1000, 1000000), 'angle': False, 'latextent': False, 'ROI': False, 'stretch': (-3,4), 'thresh': (1.0,'high'), 'text': "low"},
'ncep2-hgt-500-del2': {'area': (1000, 1000000), 'angle': False, 'latextent': False, 'ROI': False, 'stretch': (-3.5,4.5), 'thresh': (1,'high'), 'text': "depression"},
'ncep2-hgt-250-del2': {'area': (1000, 1000000), 'angle': False, 'latextent': False, 'ROI': False, 'stretch': (-6,7), 'thresh': (1,'high'), 'text': "trough"},
'ncep2-uvmag-250-del2': {'area': (1000, 1000000), 'angle': False, 'latextent': False, 'ROI': False, 'stretch': (-2.5,2.5), 'thresh': (-1,'low'), 'text': "jet"},
'ncep2-uvmag-850-del2': {'area': (1000, 1000000), 'angle': False, 'latextent': False, 'ROI': False, 'stretch': (-1.5,1.5), 'thresh': (-0.5,'low'), 'text': "jet"},
'ncep2-quvmag-700-del2': {'area': (1000, 1000000), 'angle': False, 'latextent': False, 'ROI': False, 'stretch': (-2.5,2.5), 'thresh': (-1,'low'), 'text': "moistjet"},
'ncep2-quvmag-850-del2': {'area': (1000, 1000000), 'angle': False, 'latextent': False, 'ROI': False, 'stretch': (-1.5,1.5), 'thresh': (-0.5,'low'), 'text': "moistjet"},
}
