<?xml version="1.0" encoding="utf-8"?>
<OpenGeoSysGLI>
    <name>single_fracture</name>
    <points>
        <point id="0" x="0" y="0" z="0"/>
        <point id="1" x="0" y="1" z="0"/>
        <point id="2" x="25" y="1" z="0"/>
        <point id="3" x="25" y="0" z="0"/>
        <point id="4" x="0" y="0.5" z="0" name="POINT4"/>
        <point id="5" x="25" y="0.5" z="0" name="POINT5"/>
    </points>
    <polylines>
        <polyline id="0" name="PLY_0">
            <pnt>0</pnt>
            <pnt>4</pnt>
            <pnt>5</pnt>
            <pnt>3</pnt>
            <pnt>0</pnt>
        </polyline>
        <polyline id="1" name="PLY_1">
            <pnt>1</pnt>
            <pnt>2</pnt>
            <pnt>5</pnt>
            <pnt>4</pnt>
            <pnt>1</pnt>
        </polyline>
        <polyline id="2" name="PLY_WEST">
            <pnt>0</pnt>
            <pnt>1</pnt>
        </polyline>
        <polyline id="3" name="PLY_EAST">
            <pnt>2</pnt>
            <pnt>3</pnt>
        </polyline>
        <polyline id="4" name="PLY_NORTH">
            <pnt>1</pnt>
            <pnt>2</pnt>
        </polyline>
        <polyline id="5" name="PLY_SOUTH">
            <pnt>0</pnt>
            <pnt>3</pnt>
        </polyline>
        <polyline id="6" name="PLY_FRAC">
            <pnt>4</pnt>
            <pnt>5</pnt>
        </polyline>
        <polyline id="7" name="PLY_DOMAIN">
            <pnt>0</pnt>
            <pnt>1</pnt>
            <pnt>2</pnt>
            <pnt>3</pnt>
            <pnt>0</pnt>
        </polyline>
    </polylines>
    <surfaces>
        <surface id="0" name="SURF_0">
            <element p1="0" p2="2" p3="1"/>
            <element p1="2" p2="0" p3="3"/>
        </surface>
    </surfaces>
</OpenGeoSysGLI>
