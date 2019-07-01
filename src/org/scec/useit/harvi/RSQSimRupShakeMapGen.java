package org.scec.useit.harvi;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.stat.StatUtils;
import org.opensha.commons.data.Site;
import org.opensha.commons.data.siteData.OrderedSiteDataProviderList;
import org.opensha.commons.data.siteData.SiteData;
import org.opensha.commons.data.siteData.SiteDataValue;
import org.opensha.commons.data.siteData.SiteDataValueList;
import org.opensha.commons.data.siteData.impl.CVM4i26BasinDepth;
import org.opensha.commons.data.siteData.impl.WillsMap2015;
import org.opensha.commons.data.xyz.ArbDiscrGeoDataSet;
import org.opensha.commons.data.xyz.GeoDataSet;
import org.opensha.commons.data.xyz.GriddedGeoDataSet;
import org.opensha.commons.eq.MagUtils;
import org.opensha.commons.exceptions.GMT_MapException;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationList;
import org.opensha.commons.geo.Region;
import org.opensha.commons.mapping.gmt.GMT_Map;
import org.opensha.commons.param.Parameter;
import org.opensha.commons.util.DataUtils.MinMaxAveTracker;
import org.opensha.commons.util.ExceptionUtils;
import org.opensha.commons.util.FileUtils;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.commons.util.cpt.CPTVal;
import org.opensha.sha.calc.ScenarioShakeMapCalculator;
import org.opensha.sha.earthquake.EqkRupture;
import org.opensha.sha.faultSurface.FaultTrace;
import org.opensha.sha.gui.infoTools.IMT_Info;
import org.opensha.sha.imr.AttenuationRelationship;
import org.opensha.sha.imr.attenRelImpl.NGAWest_2014_Averaged_AttenRel;
import org.opensha.sha.imr.param.IntensityMeasureParams.PGA_Param;
import org.opensha.sha.imr.param.IntensityMeasureParams.PGV_Param;
import org.opensha.sha.imr.param.IntensityMeasureParams.SA_Param;
import org.opensha.sha.mapping.GMT_MapGeneratorForShakeMaps;
import org.opensha.sha.simulators.RSQSimEvent;
import org.opensha.sha.util.ShakeMapXMLWriter;
import org.opensha.sha.util.SiteTranslator;

import com.google.common.base.Preconditions;

import scratch.UCERF3.enumTreeBranches.DeformationModels;
import scratch.UCERF3.enumTreeBranches.FaultModels;
import scratch.kevin.simulators.RSQSimCatalog;

public class RSQSimRupShakeMapGen {
	
	private final AttenuationRelationship gmpe;
	
	private final int[] eventIDs;
	private final RSQSimEvent[] events;
	private final EqkRupture[] ruptures;
	
	private final GriddedRegion reg;
	private final List<Site> sites;
	
	private ScenarioShakeMapCalculator calc = new ScenarioShakeMapCalculator();
	private OrderedSiteDataProviderList siteDataProviders;
	
//	private CPT cpt;
	
	enum RupCombineScheme {
		
		
		
		MAX_VALUE {
			@Override
			public double combine(double[] values) {
				return StatUtils.max(values);
			}
			
		},
		SUM_SQUARES {
			@Override
			public double combine(double[] values) {
				double sumSquares = 0d;
				for (double val : values)
					sumSquares += val*val;
				return Math.sqrt(sumSquares);
			}
		},
		SUM {
			@Override
			public double combine(double[] values) {
				return StatUtils.sum(values);
			}
		};
		
		public abstract double combine(double[] values);
	}
	
	public RSQSimRupShakeMapGen(RSQSimCatalog catalog, AttenuationRelationship gmpe, double rupBufferKM,
			double regionDiscretizationDegrees, int... eventIDs) throws IOException {
		this(catalog, gmpe, rupBufferKM, regionDiscretizationDegrees, catalog.loader().byIDs(eventIDs).toArray(new RSQSimEvent[0]));
	}
	
	public RSQSimRupShakeMapGen(RSQSimCatalog catalog, AttenuationRelationship gmpe, double rupBufferKM,
			double regionDiscretizationDegrees, RSQSimEvent[] events) throws IOException {
		this.gmpe = gmpe;
		Preconditions.checkArgument(events.length >= 1);
		this.events = new RSQSimEvent[events.length];
		this.ruptures = new EqkRupture[events.length];
		
		System.out.println("Building ruptures");
		this.eventIDs = new int[events.length];
		
		for (int i=0; i<events.length; i++) {
			eventIDs[i] = events[i].getID();
			ruptures[i] = catalog.getMappedSubSectRupture(events[i]);
		}
		
		// build gridded region that encompases all ruptures with the given buffer
		System.out.println("Determining region");
		this.reg = buildRegion(ruptures, rupBufferKM, regionDiscretizationDegrees);
		
		calc = new ScenarioShakeMapCalculator();
		ArrayList<SiteData<?>> providers = new ArrayList<>();
		providers.add(new WillsMap2015());
		providers.add(new CVM4i26BasinDepth(SiteData.TYPE_DEPTH_TO_2_5));
		providers.add(new CVM4i26BasinDepth(SiteData.TYPE_DEPTH_TO_1_0));
		siteDataProviders = new OrderedSiteDataProviderList(providers);
		List<SiteDataValueList<?>> siteDatas = siteDataProviders.getAllAvailableData(reg.getNodeList());
		
		// create site list
		System.out.println("Creating site list");
		sites = new ArrayList<>();
		// sets site paramters from site data above
		SiteTranslator trans = new SiteTranslator();
		LocationList nodeList = reg.getNodeList();
		for (int i=0; i<nodeList.size(); i++) {
			Location loc = nodeList.get(i);
			Site site = new Site(loc);
			
			ArrayList<SiteDataValue<?>> datas = new ArrayList<>();
			for (SiteDataValueList<?> myDatas : siteDatas)
				datas.add(myDatas.getValue(i));
			
			for (Parameter<?> param : gmpe.getSiteParams()) {
				param = (Parameter<?>) param.clone();
				trans.setParameterValue(param, datas);
				site.addParameter(param);
			}
			sites.add(site);
		}
		
//		cpt = GMT_CPT_Files.SHAKEMAP.instance();
	}
	
	private static GriddedRegion buildRegion(EqkRupture[] ruptures, double rupBufferKM, double regionDiscretizationDegrees) {
		// this wasn't working, do it the brute force way below
//		// this will be a non-rectangular region that contains all of the given ruptures, buffered by the selected rup buffer distance
//		Region combinedRegion = null;
//		for (EqkRupture rup : ruptures) {
//			FaultTrace trace = rup.getRuptureSurface().getEvenlyDiscritizedUpperEdge();
//			// remove duplicates
//			for (int i=trace.size(); --i>=1;) {
//				Location loc1 = trace.get(i);
//				Location loc2 = trace.get(i-1);
//				if (loc1.equals(loc2))
//					trace.remove(i);
//			}
//			Region bufferedRupTrace = new Region(trace, rupBufferKM);
//			if (combinedRegion == null)
//				combinedRegion = bufferedRupTrace;
//			else
//				combinedRegion = Region.union(combinedRegion, bufferedRupTrace);
//		}
//		// now make it rectangular
//		Region rectangular = new Region(new Location(combinedRegion.getMinLat(), combinedRegion.getMinLon()),
//				new Location(combinedRegion.getMaxLat(), combinedRegion.getMaxLon()));
		MinMaxAveTracker latTrack = new MinMaxAveTracker();
		MinMaxAveTracker lonTrack = new MinMaxAveTracker();
		for (EqkRupture rup : ruptures) {
			FaultTrace trace = rup.getRuptureSurface().getUpperEdge();
			for (Location loc : trace) {
				Region reg = new Region(loc, rupBufferKM);
				latTrack.addValue(reg.getMaxLat());
				latTrack.addValue(reg.getMinLat());
				lonTrack.addValue(reg.getMaxLon());
				lonTrack.addValue(reg.getMinLon());
			}
		}
		// lets make it pretty
		double minLat = getNicelyRounded(latTrack.getMin(), true, regionDiscretizationDegrees);
		double maxLat = getNicelyRounded(latTrack.getMax(), false, regionDiscretizationDegrees);
		double minLon = getNicelyRounded(lonTrack.getMin(), true, regionDiscretizationDegrees);
		double maxLon = getNicelyRounded(lonTrack.getMax(), false, regionDiscretizationDegrees);
		
		Region rectangular = new Region(new Location(minLat, minLon),
				new Location(maxLat, maxLon));
		return new GriddedRegion(rectangular, regionDiscretizationDegrees, null);
	}
	
	private static double getNicelyRounded(double val, boolean floor, double discr) {
		double numBinsAwayFromZero = Math.floor(val / discr);
		if (floor)
			numBinsAwayFromZero = Math.floor(numBinsAwayFromZero);
		else
			numBinsAwayFromZero = Math.ceil(numBinsAwayFromZero);
		return numBinsAwayFromZero * discr;
	}
	
	public void calcAll(File outputDir, RupCombineScheme... combineSchemes) throws IOException {
		// should be same order as GMT_MapGenerator.makeHazusFileSetUsingServlet
		String[] imts = { SA_Param.NAME, SA_Param.NAME, PGA_Param.NAME, PGV_Param.NAME };
		double[] periods = { 0.3, 1d, Double.NaN, Double.NaN };
		
		File[] rupDirs = new File[eventIDs.length];
		for (int i=0; i<eventIDs.length; i++) {
			rupDirs[i] = new File(outputDir, "event_"+eventIDs[i]);
			rupDirs[i].mkdir();
		}
		GriddedGeoDataSet[][] rupMaps = new GriddedGeoDataSet[eventIDs.length][imts.length];
		File[] combinedDirs;
		GriddedGeoDataSet[][] combinedMaps;
		if (eventIDs.length > 1) {
			Preconditions.checkState(combineSchemes.length > 0);
			combinedDirs = new File[combineSchemes.length];
			for (int c=0; c<combineSchemes.length; c++) {
				combinedDirs[c] = new File(outputDir, "combined_"+combineSchemes[c]);
				combinedDirs[c].mkdir();
			}
			combinedMaps = new GriddedGeoDataSet[combineSchemes.length][imts.length];
		} else {
			combinedDirs = null;
			combinedMaps = null;
		}
		
		for (int i=0; i<imts.length; i++) {
			System.out.println("Calculating for IMT "+imts[i]+" ("+periods[i]+")");
			for (int e=0; e<eventIDs.length; e++) {
				System.out.println("Calculating event "+eventIDs[e]);
				rupMaps[e][i] = calcShakeMap(ruptures[e], imts[i], periods[i]);
			}
			if (eventIDs.length > 1) {
				System.out.println("Calculating combined");
				GriddedGeoDataSet[] mapsForCombine = new GriddedGeoDataSet[eventIDs.length];
				for (int e=0; e<eventIDs.length; e++)
					mapsForCombine[e] = rupMaps[e][i];
				for (int c=0; c<combineSchemes.length; c++)
					combinedMaps[c][i] = combineShakeMaps(combineSchemes[c], mapsForCombine);
			}
		}
		
		// put it all together
		for (int e=0; e<eventIDs.length; e++)
			generateProducts(rupMaps[e], rupDirs[e], ruptures[e].getMag(), eventIDs[e]);
		if (eventIDs.length > 1) {
			double totMoment = 0;
			for (EqkRupture rup : ruptures)
				totMoment += MagUtils.magToMoment(rup.getMag());
			double combinedMag = MagUtils.momentToMag(totMoment);
			for (int c=0; c<combineSchemes.length; c++)
				generateProducts(combinedMaps[c], combinedDirs[c], combinedMag, eventIDs);
		}
	}
	
	private static final String[] extensionsToDownload = { ".shp", ".dbf", ".shx", ".png" };
	
	private void generateProducts(GriddedGeoDataSet[] datasets, File outputDir, double mag, int... ids) throws IOException {
		System.out.println("Generating maps/shape files for "+outputDir.getAbsolutePath());
		try {
			String[] prefixes = { GMT_MapGeneratorForShakeMaps.SA_03, GMT_MapGeneratorForShakeMaps.SA_10,
					GMT_MapGeneratorForShakeMaps.PGA, GMT_MapGeneratorForShakeMaps.PGV };
			GMT_Map[] maps = buildMaps(datasets, prefixes);
			String[] labels = { "0.3s SA (g)", "1s SA (g)", "PGA (g)", "PGV (cm/s)" };
			// set map file names
			for (int i=0; i<maps.length; i++) {
				maps[i].setXyzFileName(prefixes[i]+".xyz");
				maps[i].setPNGFileName(prefixes[i]+".png");
				maps[i].setCustomCptFileName(prefixes[i]+".cpt");
				maps[i].setCustomLabel(labels[i]);
			}
			String[] addresses = GMT_MapGeneratorForShakeMaps.makeHazusFileSetUsingServlet(
					maps[0], maps[1], maps[2], maps[3], "", null);
			
			String urlPrefix = addresses[0].substring(0, addresses[0].lastIndexOf("/"));
			
			for (int i = 0; i < prefixes.length; i++) {
				String prefix = prefixes[i];
				for (String ext : extensionsToDownload)
					FileUtils.downloadURL(urlPrefix+"/"+prefix+ext, new File(outputDir, prefix+ext));
				ArbDiscrGeoDataSet.writeXYZFile(maps[i].getGriddedData(), new File(outputDir, prefix+".xyz"));
			}
		} catch (GMT_MapException e) {
			ExceptionUtils.throwAsRuntimeException(e);
		}
		// now write XML file
		System.out.println("Writing XML file");
		File xmlFile = new File(outputDir, "shakemap.xml");
		String idStr;
		String name;
		if (ids.length == 1) {
			idStr = ids[0]+"";
			name = "Event "+ids[0];
		} else {
			idStr = "";
			name = "Events ";
			for (int i=0; i<eventIDs.length; i++) {
				if (i > 0) {
					name += ",";
					idStr += "_";
				}
				name += eventIDs[i];
				idStr += eventIDs[i];
			}
		}
		ShakeMapXMLWriter.writeXML(xmlFile, mag, name, idStr, datasets[2], datasets[3], datasets[0], datasets[1], null);
	}
	
	private static CPT[] getCPTs(String[] prefixes, GriddedGeoDataSet[] xyzs) {
		CPT[] cpts = new CPT[prefixes.length];
		for (int i=0; i<cpts.length; i++)
			cpts[i] = getCPT(prefixes[i], xyzs[i]);
		return cpts;
	}
	
	private static CPT getCPT(String prefix, GriddedGeoDataSet xyz) {
		List<double[]> ranges = new ArrayList<>();
		List<Color[]> colors = new ArrayList<>();
//		colors.add(toArray(new Color(255, 255, 255), new Color(255, 255, 255)));	// 0 to 1
		colors.add(toArray(new Color(255, 255, 255), new Color(191, 204, 255)));	// 1 to 2
		colors.add(toArray(new Color(191, 204, 255), new Color(160, 230, 255)));	// 2 to 3
		colors.add(toArray(new Color(160, 230, 255), new Color(128, 255, 255)));	// 3 to 4
		colors.add(toArray(new Color(128, 255, 255), new Color(122, 255, 147)));	// 4 to 5
		colors.add(toArray(new Color(122, 255, 147), new Color(255, 255, 0)));		// 5 to 6
		colors.add(toArray(new Color(255, 255, 0), new Color(255, 200, 0)));		// 6 to 7
		colors.add(toArray(new Color(255, 200, 0), new Color(255, 145, 0)));		// 7 to 8
		colors.add(toArray(new Color(255, 145, 0), new Color(255, 0, 0)));			// 8 to 9
		colors.add(toArray(new Color(255, 0, 0), new Color(200, 0, 0)));			// 9 to 10
		colors.add(toArray(new Color(200, 0, 0), new Color(128, 0, 0)));			// 10 to 13
		switch (prefix) {
		case GMT_MapGeneratorForShakeMaps.PGA:
//			ranges.add(toArray(0d, 0d));
			ranges.add(toArray(0d, 0.17));
			ranges.add(toArray(0.17, logMean(0.17, 1.4)));
			ranges.add(toArray(logMean(0.17, 1.4), 1.4));
			ranges.add(toArray(1.4, 3.9));
			ranges.add(toArray(3.9, 9.2));
			ranges.add(toArray(9.2, 18));
			ranges.add(toArray(18, 34));
			ranges.add(toArray(34, 65));
			ranges.add(toArray(65, 124));
			ranges.add(toArray(124, 250));
			scale(ranges, 1d/100d); // convert to g
			break;
		case GMT_MapGeneratorForShakeMaps.PGV:
//			ranges.add(toArray(0d, 0d));
			ranges.add(toArray(0d, 0.1));
			ranges.add(toArray(0.1, logMean(0.1, 1.1)));
			ranges.add(toArray(logMean(0.1, 1.1), 1.1));
			ranges.add(toArray(1.1, 3.4));
			ranges.add(toArray(3.4, 8.1));
			ranges.add(toArray(8.1, 16));
			ranges.add(toArray(16, 31));
			ranges.add(toArray(31, 60));
			ranges.add(toArray(60, 116));
			ranges.add(toArray(116, 250));
			break;

		default:
			break;
		}
		if (ranges.isEmpty()) {
			// log scale between 1e-3 and max value
			double min = 1e-3;
			double max = xyz.getMaxZ();
			double logMin = Math.log(min);
			double logMax = Math.log(max);
			double logPer = (logMax - logMin)/(colors.size());
			for (int i=0; i<colors.size(); i++) {
				ranges.add(toArray(Math.exp(logMin + logPer*i), Math.exp(logMin + logPer*(i+1))));
			}
		}
		Preconditions.checkState(ranges.size() == colors.size());
		
		CPT cpt = new CPT();
		
		for (int i=0; i<colors.size(); i++) {
			double[] range = ranges.get(i);
			Color[] c = colors.get(i);
			CPTVal val = new CPTVal((float)range[0], c[0], (float)range[1], c[1]);
			cpt.add(val);
		}
		cpt.setBelowMinColor(cpt.getMinColor());
		cpt.setAboveMaxColor(cpt.getMaxColor());
		
		return cpt;
	}
	
	private static void scale(List<double[]> ranges, double scalar) {
		for (int i=0; i<ranges.size(); i++) {
			double[] vals = ranges.get(i);
			for (int j=0; j<vals.length; j++)
				vals[j] *= scalar;
		}
	}
	
	private static double logMean(double val1, double val2) {
		return Math.exp(0.5*Math.log(val1) + 0.5*Math.log(val2));
	}
	
	private static Color[] toArray(Color... colors) {
		return colors;
	}
	private static double[] toArray(double... vals) {
		return vals;
	}
	
	private GMT_Map[] buildMaps(GriddedGeoDataSet[] xyzs, String[] prefixes) {
		GMT_Map[] maps = new GMT_Map[xyzs.length];
		for (int i=0; i<xyzs.length; i++)
			maps[i] = buildMap(xyzs[i], prefixes[i]);
		return maps;
	}
	
	private GMT_Map buildMap(GriddedGeoDataSet xyz, String prefix) {
////		CPT cpt = this.cpt.rescale(0d, Math.exp(xyz.getMaxZ()));
//		CPT cpt = this.cpt.rescale(0d, xyz.getMaxZ());
		CPT cpt = getCPT(prefix, xyz);
		GMT_Map map = new GMT_Map(reg, xyz, reg.getLatSpacing(), cpt);
//		map.setCustomScaleMin((double)cpt.getMinValue());
//		map.setCustomScaleMax((double)cpt.getMaxValue());
		map.setCustomCptFileName(prefix+".cpt");
		map.setRescaleCPT(false);
		map.setLogPlot(false);
		map.setUseGMTSmoothing(true);
		return map;
	}
	
	public GriddedGeoDataSet calcShakeMap(String imt, double period, RupCombineScheme combinationScheme) {
		if (ruptures.length > 1)
			Preconditions.checkNotNull(combinationScheme);
		
		GriddedGeoDataSet[] maps = new GriddedGeoDataSet[ruptures.length];
		
		for (int r=0; r<ruptures.length; r++) {
			EqkRupture rup = ruptures[r];
			maps[r] = calcShakeMap(rup, imt, period);
		}
		
		return combineShakeMaps(combinationScheme, maps);
	}
	
	public GriddedGeoDataSet combineShakeMaps(RupCombineScheme combinationScheme, GriddedGeoDataSet... maps) {
		if (maps.length == 1)
			return maps[0];
		
		// more than one, must combine
		GriddedGeoDataSet combined = new GriddedGeoDataSet(reg, maps[0].isLatitudeX());
		
		for (int i=0; i<maps[0].size(); i++) {
			double[] vals = new double[maps.length];
			for (int j=0; j<maps.length; j++) {
				GriddedGeoDataSet map = maps[j];
				vals[j] = map.get(i);
			}
			combined.set(maps[0].getLocation(i), combinationScheme.combine(vals));
		}
		
		return combined;
	}
	
	public GriddedGeoDataSet calcShakeMap(EqkRupture rup, String imt, double period) {
		gmpe.setIntensityMeasure(imt);
		if (imt.equals(SA_Param.NAME))
			SA_Param.setPeriodInSA_Param(gmpe.getIntensityMeasure(), period);
		List<AttenuationRelationship> selectedAttenRels = new ArrayList<>();
		List<Double> attenRelWts = new ArrayList<>();
		selectedAttenRels.add(gmpe);
		attenRelWts.add(1d);
		GeoDataSet xyz = calc.getScenarioShakeMapData(selectedAttenRels, attenRelWts, sites, rup, false, 0.5);
		GriddedGeoDataSet gridXYZ = new GriddedGeoDataSet(reg, xyz.isLatitudeX());
		for (int i=0; i<xyz.size(); i++)
			gridXYZ.set(xyz.getLocation(i), xyz.get(i));
		if (IMT_Info.isIMT_LogNormalDist(imt))
			gridXYZ.exp();
		System.out.println("Done calculating, data range: "+gridXYZ.getMinZ()+" to "+gridXYZ.getMaxZ());
		return gridXYZ;
	}
	
	public static void main(String[] args) throws IOException {
		// directory containing all input files, including the .flt file and the e/p/t/dList files
		// change this to your directory
		File catalogDir = new File("/data/kevin/simulators/catalogs/bruce/rundir2585");
		RSQSimCatalog catalog = new RSQSimCatalog(catalogDir, "Catalog", FaultModels.FM3_1, DeformationModels.GEOLOGIC);
		// output directory will files will be saved. change this to a sensible place for you
//		File outputDir = new File("C:\\Users\\scecintern2016\\Desktop\\Catalog\\output");
		File outputDir = new File("/tmp/shakemap");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		// event ID(s) to calculate shakemaps for. if more than one, then they will be combined using the combining schemes given below
		int[] eventIDs = {1056770};
//		int[] eventIDs = { 956271, 956278};
		// scheme(s) for combining multiple events
		RupCombineScheme[] combiningSchemes = { RupCombineScheme.MAX_VALUE, RupCombineScheme.SUM_SQUARES, RupCombineScheme.SUM };
		// buffer around ruptures in kilometers
		double rupBufferKM = 150;
		// Map discretization in degrees. Set low (e.g. 0.02 to 0.1) for nice looking high resolution maps that take long to calculate.
		// Set high (e.g. 0.1 to 0.5) for ugly maps that calculate quickly, mostly for testing purposes.
		// It's an exponential scale, so 0.1 degrees will take 4 times longer than 0.2 degrees,
		// and 0.05 degrees 4 times longer than 0.1 degrees, etc...
		// speed also depends on the size of the region, which is determined by the rup buffer size and the size of the rupture itself
		double regionDiscretizationDegrees = 0.05;
		
		System.out.println("Loading catalog from "+catalogDir.getAbsolutePath());
		List<RSQSimEvent> events = catalog.loader().byIDs(eventIDs);
		System.out.println("Loaded "+events.size()+" events");
		
		AttenuationRelationship gmpe = new NGAWest_2014_Averaged_AttenRel(null, false);
		gmpe.setParamDefaults();
		
		gmpe.setIntensityMeasure(PGA_Param.NAME);
		
		RSQSimRupShakeMapGen mapGen = new RSQSimRupShakeMapGen(catalog, gmpe, rupBufferKM,
				regionDiscretizationDegrees, events.toArray(new RSQSimEvent[0]));
		
		mapGen.calcAll(outputDir, combiningSchemes);
	}

}
