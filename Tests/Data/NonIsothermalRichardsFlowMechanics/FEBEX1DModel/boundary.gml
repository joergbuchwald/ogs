<?xml version="1.0" encoding="ISO-8859-1"?>
<?xml-stylesheet type="text/xsl" href="OpenGeoSysGLI.xsl"?>

<OpenGeoSysGLI xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xmlns:ogs="http://www.opengeosys.org">
    <name>febex1D_geo</name>
    <points>
        <point id="0" x="3.42" y="0.0" z="0"/>
        <point id="1" x="3.42" y="0.02" z="0"/>
        <point id="2" x="0.45" y="0.0" z="0"/>
        <point id="3" x="0.45" y="0.02" z="0"/>
    </points>

    <polylines>
        <polyline id="1" name="right">
            <pnt>0</pnt>
            <pnt>1</pnt>
        </polyline>
        <polyline id="2" name="bottom">
            <pnt>2</pnt>
            <pnt>0</pnt>
        </polyline>
        <polyline id="3" name="top">
            <pnt>3</pnt>
            <pnt>1</pnt>
        </polyline>
        <polyline id="4" name="heater_surface">
            <pnt>2</pnt>
            <pnt>3</pnt>
        </polyline>
    </polylines>
</OpenGeoSysGLI>
