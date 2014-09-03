/* ========================================================================
    File: @(#)gdal_nodedensity.cpp
   ------------------------------------------------------------------------
    plot OSM node/way density with gdal
    based on gdal_rasterize by Frank Warmerdam and nodedensity by Jochen Topf
   ------------------------------------------------------------------------

 ******************************************************************************
 * Copyright (c) 2005, Frank Warmerdam <warmerdam@pobox.com>
 * Copyright (c) 2013, Jochen Topf <jochen@topf.org>
 * Copyright (c) 2014, Christoph Hormann <chris_hormann@gmx.de>
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the "Software"),
 * to deal in the Software without restriction, including without limitation
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,
 * and/or sell copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
 * THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 * DEALINGS IN THE SOFTWARE.
 ****************************************************************************/

#include "gdal.h"
#include "gdal_alg.h"
#include "cpl_conv.h"
#include "ogr_api.h"
#include "ogr_srs_api.h"
#include "cpl_string.h"
#include "commonutils.h"
#include <vector>

#define OSMIUM_WITH_PBF_INPUT
#define OSMIUM_WITH_XML_INPUT

#include <osmium.hpp>
#include <osmium/storage/byid/sparse_table.hpp>
#include <osmium/storage/byid/mmap_file.hpp>
#include <osmium/handler/coordinates_for_ways.hpp>
#include <osmium/geometry/point.hpp>
#include <osmium/geometry/ogr.hpp>

// data type for the density counters
typedef unsigned int CounterType;
// corresponding GDAL data type
#define CounterTypeGDT GDT_UInt32

/************************************************************************/
/*                            ArgIsNumeric()                            */
/************************************************************************/

static int ArgIsNumeric( const char *pszArg )
{
    return CPLGetValueType(pszArg) != CPL_VALUE_STRING;
}

/************************************************************************/
/*                               Usage()                                */
/************************************************************************/

static void Usage()
{
    printf( 
        "Usage: gdal_nodedensity [-b band]*\n"
        "       [-of format] [-a_srs srs_def] [-co \"NAME=VALUE\"]*\n"
        "       [-interval dist] [-a_nodata value] [-init value]*\n"
        "       [-te xmin ymin xmax ymax] [-tr xres yres] [-tap] [-ts width height]\n"
        "       [-ot {Byte/Int16/UInt16/UInt32/Int32/Float32/Float64/\n"
        "             CInt16/CInt32/CFloat32/CFloat64}] [-q]\n"
        "       <src_datasource> <dst_filename>\n" );
    exit( 1 );
}


/* ================================================== */

class NodeDensityHandler : public Osmium::Handler::Base {

	std::vector<CounterType*> m_Data;
	int m_xsize;
	int m_ysize;
	OGRCoordinateTransformationH m_hCT;
	double *m_GeoTransform;
	double m_InvGeoTransform[6];
	size_t m_count_all;
	size_t m_count;

public:

	NodeDensityHandler(std::vector<CounterType*> pData, int xSize, int ySize, OGRCoordinateTransformationH hCT, double *GeoTransform)
		: Base(),
			m_GeoTransform(GeoTransform),
			m_hCT(hCT),
			m_xsize(xSize),
			m_ysize(ySize),
			m_Data(pData) { 

		GDALInvGeoTransform(m_GeoTransform, m_InvGeoTransform);
		m_count = 0;
		m_count_all = 0;
	}

	~NodeDensityHandler() { }

	void node(const shared_ptr<Osmium::OSM::Node const>& node) {

		if (!node->position().defined()) return;

		double x = node->position().lon();
		double y = node->position().lat();

		m_count_all++;

		if (!OCTTransform( m_hCT, 1, &x, &y, NULL ))
			return;

		int iPixel, iLine;

		iPixel = (int) floor(m_InvGeoTransform[0] 
												 + m_InvGeoTransform[1] * x
												 + m_InvGeoTransform[2] * y);
		iLine = (int) floor(m_InvGeoTransform[3] 
												+ m_InvGeoTransform[4] * x
												+ m_InvGeoTransform[5] * y);

		if (iPixel < 0) return;
		if (iLine < 0) return;
		if (iPixel >= m_xsize) return;
		if (iLine >= m_ysize) return;

		for( unsigned int iBand = 0; iBand < m_Data.size(); iBand++ )
		{
			m_Data[iBand][iLine*m_xsize + iPixel]++;
		}

		m_count++;

		if ((m_count % 10000) == 0)
		{
			fprintf(stderr, "         \r %d/%d nodes (inside/total)...", m_count, m_count_all);
		}
	}

	void after_nodes() {
		fprintf(stderr, "\rprocessed %d/%d nodes.            \n", m_count, m_count_all);
		throw Osmium::Handler::StopReading();
	}

}; // class NodeDensityHandler

typedef Osmium::Storage::ById::SparseTable<Osmium::OSM::Position> storage_sparsetable_t;
typedef Osmium::Storage::ById::MmapFile<Osmium::OSM::Position> storage_mmap_t;
typedef Osmium::Handler::CoordinatesForWays<storage_sparsetable_t, storage_mmap_t> cfw_handler_t;

class WayDensityHandler : public Osmium::Handler::Base {

	storage_sparsetable_t store_pos;
	storage_mmap_t store_neg;
	cfw_handler_t* handler_cfw;

	std::vector<CounterType*> m_Data;
	int m_xsize;
	int m_ysize;
	OGRCoordinateTransformationH m_hCT;
	double *m_GeoTransform;
	double m_InvGeoTransform[6];
	double m_interval;
	size_t m_count_all;
	size_t m_count;
	size_t m_count_node;
	size_t m_px_count;
	double dsum;

public:

	WayDensityHandler(std::vector<CounterType*> pData, int xSize, int ySize, OGRCoordinateTransformationH hCT, double *GeoTransform, double Interval)
		: Base(),
			m_interval(Interval),
			m_GeoTransform(GeoTransform),
			m_hCT(hCT),
			m_xsize(xSize),
			m_ysize(ySize),
			m_Data(pData) { 

		handler_cfw = new cfw_handler_t(store_pos, store_neg);

		GDALInvGeoTransform(m_GeoTransform, m_InvGeoTransform);
		m_count = 0;
		m_count_node = 0;
		m_count_all = 0;
		m_px_count = 0;
		dsum = 0.0;
	}

	~WayDensityHandler() {
		delete handler_cfw;
	}

	void init(Osmium::OSM::Meta& meta) {
		handler_cfw->init(meta);
	}

	void node(const shared_ptr<Osmium::OSM::Node const>& node) {
		handler_cfw->node(node);
		m_count_node++;
		if ((m_count_node % 10000) == 0)
		{
			fprintf(stderr, "         \r %d nodes...", m_count_node);
		}
	}

	void after_nodes() {
		fprintf(stderr, "         \rprocessed %d nodes.\n", m_count_node);
		std::cerr << "\nMemory used for node coordinates storage (approximate):\n  for positive IDs: "
							<< store_pos.used_memory() / (1024 * 1024)
							<< " MiB\n  for negative IDs: "
							<< store_neg.used_memory() / (1024 * 1024)
							<< " MiB\n";

		handler_cfw->after_nodes();
	}

	void before_ways() const {
		if (m_count_node == 0)
			throw Osmium::Handler::StopReading();
	}

	void way(const shared_ptr<Osmium::OSM::Way>& way) {
		handler_cfw->way(way);

		Osmium::OSM::WayNodeList::const_iterator iter = way->nodes().begin();
		Osmium::OSM::WayNodeList::const_iterator iter2 = iter;
		iter2++;
		for (; iter2 != way->nodes().end(); ++iter, ++iter2)
		{
			if (!iter->position().defined()) continue;
			if (!iter2->position().defined()) continue;

			double x = 0.5*(iter->lon()+iter2->lon());
			double y = 0.5*(iter->lat()+iter2->lat());

			m_count_all++;

			if (!OCTTransform( m_hCT, 1, &x, &y, NULL ))
				return;

			int iPixel, iLine;

			iPixel = (int) floor(m_InvGeoTransform[0] 
													 + m_InvGeoTransform[1] * x
													 + m_InvGeoTransform[2] * y);
			iLine = (int) floor(m_InvGeoTransform[3] 
													+ m_InvGeoTransform[4] * x
													+ m_InvGeoTransform[5] * y);

			if (iPixel < 0) continue;
			if (iLine < 0) continue;
			if (iPixel >= m_xsize) continue;
			if (iLine >= m_ysize) continue;

			double x2 = iter->lon();
			double y2 = iter->lat();

			if (!OCTTransform( m_hCT, 1, &x2, &y2, NULL ))
				return;

			double dx = 2.0*(x2-x);
			double dy = 2.0*(y2-y);

			double d = std::sqrt(dx*dx + dy*dy)/m_interval;

			for( unsigned int iBand = 0; iBand < m_Data.size(); iBand++ )
			{
				int di = (int)floor(d+0.5);
				if (di < 1) di = 1;
				dsum += di;
				m_px_count++;
				m_Data[iBand][iLine*m_xsize + iPixel] += di;
			}

			m_count++;

			if ((m_count % 10000) == 0)
			{
				fprintf(stderr, "         \r %d/%d segments (inside/total)...", m_count, m_count_all);
			}
		}
	}

	void after_ways() {
		fprintf(stderr, "\rprocessed %d/%d nodes (average %.2f segments).      \n", m_count, m_count_all, dsum/m_px_count);
		throw Osmium::Handler::StopReading();
	}
};

/* ================================================== */

static void ProcessOSM(const char *osm_fnm, GDALDatasetH hDstDS, std::vector<int> anBandList, const double Interval)
{
	std::vector<CounterType*> pData;
	int nXSize, nYSize;

	for( unsigned int iBand = 0; iBand < anBandList.size(); iBand++ )
	{
		GDALRasterBandH hBand = GDALGetRasterBand( hDstDS, iBand+1 );

		nXSize = GDALGetRasterBandXSize(hBand);
		nYSize = GDALGetRasterBandYSize(hBand);

		pData.push_back((CounterType *) CPLMalloc(sizeof(CounterType)*nXSize*nYSize));
	}

	OGRSpatialReferenceH hSrcSRS = NULL;
	hSrcSRS = OSRNewSpatialReference(NULL);
	OSRImportFromEPSG( hSrcSRS, 4326 );
	OGRSpatialReferenceH  hDstSRS = NULL;

	if( GDALGetProjectionRef( hDstDS ) != NULL )
	{
		char *pszProjection;
    
		pszProjection = (char *) GDALGetProjectionRef( hDstDS );
    
		hDstSRS = OSRNewSpatialReference(NULL);
		if( OSRImportFromWkt( hDstSRS, &pszProjection ) != CE_None )
		{
			OSRDestroySpatialReference(hDstSRS);
			hDstSRS = NULL;
		}
	}
	else
	{
		hDstSRS = OSRNewSpatialReference(NULL);
		OSRImportFromEPSG( hDstSRS, 4326 );
	}
 
	OGRCoordinateTransformationH hCT = OCTNewCoordinateTransformation(hSrcSRS, hDstSRS);

	if (hCT == NULL)
	{
		fprintf(stderr, "Failed to create coordinate system transform.\n");
		return;
	}

	double adfGeoTransform[6];
	if( GDALGetGeoTransform( hDstDS, adfGeoTransform ) != CE_None )
	{
		fprintf(stderr, "Failed to get geotransform from raster.\n");
		return;
	}

	if (Interval <= 0.0)
	{
		Osmium::OSMFile infile(osm_fnm);
		NodeDensityHandler handler(pData, nXSize, nYSize, hCT, adfGeoTransform);
		Osmium::Input::read(infile, handler);
	}
	else
	{
		Osmium::OSMFile infile(osm_fnm);
		WayDensityHandler handler(pData, nXSize, nYSize, hCT, adfGeoTransform, Interval);
		Osmium::Input::read(infile, handler);
	}

	OCTDestroyCoordinateTransformation(hCT);

	for( unsigned int iBand = 0; iBand < anBandList.size(); iBand++ )
	{
		GDALRasterBandH hBand = GDALGetRasterBand( hDstDS, iBand+1 );
		if (GDALGetRasterDataType(hBand) == GDT_Byte)
		{
			for (int i = 0; i<nXSize*nYSize; i++)
			{
				if (pData[iBand][i] > 255) pData[iBand][i] = 255;
			}
		}
		else if (GDALGetRasterDataType(hBand) == GDT_UInt16)
		{
			for (int i = 0; i<nXSize*nYSize; i++)
			{
				if (pData[iBand][i] > 65535) pData[iBand][i] = 65535;
			}
		}

		GDALRasterIO(hBand, GF_Write, 0, 0, nXSize, nYSize, pData[iBand], nXSize, nYSize, CounterTypeGDT, 0, 0);
		CPLFree( pData[iBand] );
	}
}

/************************************************************************/
/*                  CreateOutputDataset()                               */
/************************************************************************/

static
GDALDatasetH CreateOutputDataset(
                                 OGRSpatialReferenceH hSRS,
                                 int bGotBounds, OGREnvelope sEnvelop,
                                 GDALDriverH hDriver, const char* pszDstFilename,
                                 int nXSize, int nYSize, double dfXRes, double dfYRes,
                                 int bTargetAlignedPixels,
                                 int nBandCount, GDALDataType eOutputType,
                                 char** papszCreateOptions, std::vector<double> adfInitVals,
                                 int bNoDataSet, double dfNoData)
{
    int bFirstLayer = TRUE;
    char* pszWKT = NULL;
    GDALDatasetH hDstDS = NULL;
    unsigned int i;

    if (dfXRes == 0 && dfYRes == 0)
    {
        dfXRes = (sEnvelop.MaxX - sEnvelop.MinX) / nXSize;
        dfYRes = (sEnvelop.MaxY - sEnvelop.MinY) / nYSize;
    }
    else if (bTargetAlignedPixels && dfXRes != 0 && dfYRes != 0)
    {
        sEnvelop.MinX = floor(sEnvelop.MinX / dfXRes) * dfXRes;
        sEnvelop.MaxX = ceil(sEnvelop.MaxX / dfXRes) * dfXRes;
        sEnvelop.MinY = floor(sEnvelop.MinY / dfYRes) * dfYRes;
        sEnvelop.MaxY = ceil(sEnvelop.MaxY / dfYRes) * dfYRes;
    }

    double adfProjection[6];
    adfProjection[0] = sEnvelop.MinX;
    adfProjection[1] = dfXRes;
    adfProjection[2] = 0;
    adfProjection[3] = sEnvelop.MaxY;
    adfProjection[4] = 0;
    adfProjection[5] = -dfYRes;

    if (nXSize == 0 && nYSize == 0)
    {
        nXSize = (int)(0.5 + (sEnvelop.MaxX - sEnvelop.MinX) / dfXRes);
        nYSize = (int)(0.5 + (sEnvelop.MaxY - sEnvelop.MinY) / dfYRes);
    }

    hDstDS = GDALCreate(hDriver, pszDstFilename, nXSize, nYSize,
                        nBandCount, eOutputType, papszCreateOptions);
    if (hDstDS == NULL)
    {
        fprintf(stderr, "Cannot create %s\n", pszDstFilename);
        exit(2);
    }

    GDALSetGeoTransform(hDstDS, adfProjection);

    if (hSRS)
        OSRExportToWkt(hSRS, &pszWKT);
    if (pszWKT)
        GDALSetProjection(hDstDS, pszWKT);
    CPLFree(pszWKT);

    int iBand;

    if (bNoDataSet)
    {
        for(iBand = 0; iBand < nBandCount; iBand++)
        {
            GDALRasterBandH hBand = GDALGetRasterBand(hDstDS, iBand + 1);
            GDALSetRasterNoDataValue(hBand, dfNoData);
        }
    }

    if (adfInitVals.size() != 0)
    {
        for(iBand = 0; iBand < MIN(nBandCount,(int)adfInitVals.size()); iBand++)
        {
            GDALRasterBandH hBand = GDALGetRasterBand(hDstDS, iBand + 1);
            GDALFillRaster(hBand, adfInitVals[iBand], 0);
        }
    }

    return hDstDS;
}

/************************************************************************/
/*                                main()                                */
/************************************************************************/

int main( int argc, char ** argv )

{
    int i = FALSE;
    const char *pszSrcFilename = NULL;
    const char *pszDstFilename = NULL;
    std::vector<int> anBandList;
    char **papszRasterizeOptions = NULL;
    double dfXRes = 0, dfYRes = 0;
    int bCreateOutput = FALSE;
    const char* pszFormat = "GTiff";
    int bFormatExplicitelySet = FALSE;
    char **papszCreateOptions = NULL;
    GDALDriverH hDriver = NULL;
    GDALDataType eOutputType = GDT_Byte;
    std::vector<double> adfInitVals;
    int bNoDataSet = FALSE;
    double dfNoData = 0;
    OGREnvelope sEnvelop;
    int bGotBounds = FALSE;
    int nXSize = 0, nYSize = 0;
    int bQuiet = FALSE;
    GDALProgressFunc pfnProgress = GDALTermProgress;
    OGRSpatialReferenceH hSRS = NULL;
    int bTargetAlignedPixels = FALSE;
		double Interval = -1.0;
    

    /* Check that we are running against at least GDAL 1.4 */
    /* Note to developers : if we use newer API, please change the requirement */
    if (atoi(GDALVersionInfo("VERSION_NUM")) < 1400)
    {
        fprintf(stderr, "At least, GDAL >= 1.4.0 is required for this version of %s, "
                "which was compiled against GDAL %s\n", argv[0], GDAL_RELEASE_NAME);
        exit(1);
    }

    GDALAllRegister();
    OGRRegisterAll();

    argc = GDALGeneralCmdLineProcessor( argc, &argv, 0 );
    if( argc < 1 )
        exit( -argc );

/* -------------------------------------------------------------------- */
/*      Parse arguments.                                                */
/* -------------------------------------------------------------------- */
    for( i = 1; i < argc; i++ )
    {
        if( EQUAL(argv[i], "--utility_version") )
        {
            printf("%s was compiled against GDAL %s and is running against GDAL %s\n",
                   argv[0], GDAL_RELEASE_NAME, GDALVersionInfo("RELEASE_NAME"));
            return 0;
        }
        else if( EQUAL(argv[i],"-q") || EQUAL(argv[i],"-quiet") )
        {
            bQuiet = TRUE;
            pfnProgress = GDALDummyProgress;
        }
        else if( EQUAL(argv[i],"-b") && i < argc-1 )
        {
            if (strchr(argv[i+1], ' '))
            {
                char** papszTokens = CSLTokenizeString( argv[i+1] );
                char** papszIter = papszTokens;
                while(papszIter && *papszIter)
                {
                    anBandList.push_back(atoi(*papszIter));
                    papszIter ++;
                }
                CSLDestroy(papszTokens);
                i += 1;
            }
            else
            {
                while(i < argc-1 && ArgIsNumeric(argv[i+1]))
                {
                    anBandList.push_back(atoi(argv[i+1]));
                    i += 1;
                }
            }
        }
        else if( EQUAL(argv[i],"-of") && i < argc-1 )
        {
            pszFormat = argv[++i];
            bFormatExplicitelySet = TRUE;
            bCreateOutput = TRUE;
        }
        else if( EQUAL(argv[i],"-a_nodata") && i < argc - 1 )
        {
            dfNoData = atof(argv[i+1]);
            bNoDataSet = TRUE;
            i += 1;
            bCreateOutput = TRUE;
        }
        else if( EQUAL(argv[i],"-a_srs") && i < argc-1 )
        {
            hSRS = OSRNewSpatialReference( NULL );

            if( OSRSetFromUserInput(hSRS, argv[i+1]) != OGRERR_NONE )
            {
                fprintf( stderr, "Failed to process SRS definition: %s\n", 
                         argv[i+1] );
                exit( 1 );
            }

            i++;
            bCreateOutput = TRUE;
        }   
        else if( EQUAL(argv[i],"-interval") && i < argc - 1 )
        {
            Interval = atof(argv[++i]);
        }
        else if( EQUAL(argv[i],"-te") && i < argc - 4 )
        {
            sEnvelop.MinX = atof(argv[++i]);
            sEnvelop.MinY = atof(argv[++i]);
            sEnvelop.MaxX = atof(argv[++i]);
            sEnvelop.MaxY = atof(argv[++i]);
            bGotBounds = TRUE;
            bCreateOutput = TRUE;
        }
        else if( EQUAL(argv[i],"-a_ullr") && i < argc - 4 )
        {
            sEnvelop.MinX = atof(argv[++i]);
            sEnvelop.MaxY = atof(argv[++i]);
            sEnvelop.MaxX = atof(argv[++i]);
            sEnvelop.MinY = atof(argv[++i]);
            bGotBounds = TRUE;
            bCreateOutput = TRUE;
        }
        else if( EQUAL(argv[i],"-co") && i < argc-1 )
        {
            papszCreateOptions = CSLAddString( papszCreateOptions, argv[++i] );
            bCreateOutput = TRUE;
        }
        else if( EQUAL(argv[i],"-ot") && i < argc-1 )
        {
            int	iType;
            
            for( iType = 1; iType < GDT_TypeCount; iType++ )
            {
                if( GDALGetDataTypeName((GDALDataType)iType) != NULL
                    && EQUAL(GDALGetDataTypeName((GDALDataType)iType),
                             argv[i+1]) )
                {
                    eOutputType = (GDALDataType) iType;
                }
            }

            if( eOutputType == GDT_Unknown )
            {
                printf( "Unknown output pixel type: %s\n", argv[i+1] );
                Usage();
            }
            i++;
            bCreateOutput = TRUE;
        }
        else if( (EQUAL(argv[i],"-ts") || EQUAL(argv[i],"-outsize")) && i < argc-2 )
        {
            nXSize = atoi(argv[++i]);
            nYSize = atoi(argv[++i]);
            if (nXSize <= 0 || nYSize <= 0)
            {
                printf( "Wrong value for -outsize parameters\n");
                Usage();
            }
            bCreateOutput = TRUE;
        }
        else if( EQUAL(argv[i],"-tr") && i < argc-2 )
        {
            dfXRes = atof(argv[++i]);
            dfYRes = fabs(atof(argv[++i]));
            if( dfXRes == 0 || dfYRes == 0 )
            {
                printf( "Wrong value for -tr parameters\n");
                Usage();
            }
            bCreateOutput = TRUE;
        }
        else if( EQUAL(argv[i],"-tap") )
        {
            bTargetAlignedPixels = TRUE;
            bCreateOutput = TRUE;
        }
        else if( pszSrcFilename == NULL )
        {
            pszSrcFilename = argv[i];
        }
        else if( pszDstFilename == NULL )
        {
            pszDstFilename = argv[i];
        }
        else
            Usage();
    }

    if( pszSrcFilename == NULL || pszDstFilename == NULL )
    {
        fprintf( stderr, "Missing source or destination.\n\n" );
        Usage();
    }

    if( bCreateOutput )
    {
        if( dfXRes == 0 && dfYRes == 0 && nXSize == 0 && nYSize == 0 )
        {
            fprintf( stderr, "'-tr xres yes' or '-ts xsize ysize' is required.\n\n" );
            Usage();
        }
    
        if (bTargetAlignedPixels && dfXRes == 0 && dfYRes == 0)
        {
            fprintf( stderr, "-tap option cannot be used without using -tr\n");
            Usage();
        }

        if( anBandList.size() != 0 )
        {
            fprintf( stderr, "-b option cannot be used when creating a GDAL dataset.\n\n" );
            Usage();
        }

        int nBandCount = 1;

        int i;
        for(i=1;i<=nBandCount;i++)
            anBandList.push_back( i );
    }
    else
    {
        if( anBandList.size() == 0 )
            anBandList.push_back( 1 );
    }


/* -------------------------------------------------------------------- */
/*      Open target raster file.  Eventually we will add optional       */
/*      creation.                                                       */
/* -------------------------------------------------------------------- */
    GDALDatasetH hDstDS = NULL;

    if (bCreateOutput)
    {
/* -------------------------------------------------------------------- */
/*      Find the output driver.                                         */
/* -------------------------------------------------------------------- */
        hDriver = GDALGetDriverByName( pszFormat );
        if( hDriver == NULL 
            || GDALGetMetadataItem( hDriver, GDAL_DCAP_CREATE, NULL ) == NULL )
        {
            int	iDr;

            printf( "Output driver `%s' not recognised or does not support\n", 
                    pszFormat );
            printf( "direct output file creation.  The following format drivers are configured\n"
                    "and support direct output:\n" );

            for( iDr = 0; iDr < GDALGetDriverCount(); iDr++ )
            {
                GDALDriverH hDriver = GDALGetDriver(iDr);

                if( GDALGetMetadataItem( hDriver, GDAL_DCAP_CREATE, NULL) != NULL )
                {
                    printf( "  %s: %s\n",
                            GDALGetDriverShortName( hDriver  ),
                            GDALGetDriverLongName( hDriver ) );
                }
            }
            printf( "\n" );
            exit( 1 );
        }

        if (!bQuiet && !bFormatExplicitelySet)
            CheckExtensionConsistency(pszDstFilename, pszFormat);
    }
    else
    {
        hDstDS = GDALOpen( pszDstFilename, GA_Update );
        if( hDstDS == NULL )
            exit( 2 );
    }

/* -------------------------------------------------------------------- */
/*      Create output file if necessary.                                */
/* -------------------------------------------------------------------- */

    if (bCreateOutput && hDstDS == NULL)
    {
        hDstDS = CreateOutputDataset(hSRS,
                                bGotBounds, sEnvelop,
                                hDriver, pszDstFilename,
                                nXSize, nYSize, dfXRes, dfYRes,
                                bTargetAlignedPixels,
                                anBandList.size(), eOutputType,
                                papszCreateOptions, adfInitVals,
                                bNoDataSet, dfNoData);
    }

/* -------------------------------------------------------------------- */
/*      Process osm Data                                                */
/* -------------------------------------------------------------------- */

		ProcessOSM(pszSrcFilename, hDstDS, anBandList, Interval);

/* -------------------------------------------------------------------- */
/*      Cleanup                                                         */
/* -------------------------------------------------------------------- */

    GDALClose( hDstDS );

    OSRDestroySpatialReference(hSRS);

    CSLDestroy( argv );
    CSLDestroy( papszRasterizeOptions );
    CSLDestroy( papszCreateOptions );
    
    GDALDestroyDriverManager();
    OGRCleanupAll();

    google::protobuf::ShutdownProtobufLibrary();

    return 0;
}
