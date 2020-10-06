<?php


namespace App\Classes;


use InvalidArgumentException;
use UnexpectedValueException;

/**
 * This class was created based on the source code of TurfJS
 * https://npmcdn.com/@turf/turf@5.1.6/turf.js
 * The purpose mainly focuses on calculate how many Qude/tile (each square meters) in a polygon
 * I did not get every functions from TurfJS.
 **/
class TurfJS {
	public $earthRadius = 6371008.8;

	public $factors = [];

	public function __construct() {
		$this->factors = [
			"meters" => $this->earthRadius,
			"metres" => $this->earthRadius,
			"millimeters" => $this->earthRadius * 1000,
			"millimetres" => $this->earthRadius * 1000,
			"centimeters" => $this->earthRadius * 100,
			"centimetres" => $this->earthRadius * 100,
			"kilometers" => $this->earthRadius / 1000,
			"kilometres" => $this->earthRadius / 1000,
			"miles" => $this->earthRadius / 1609.344,
			"nauticalmiles" => $this->earthRadius / 1852,
			"inches" => $this->earthRadius * 39.370,
			"yards" => $this->earthRadius / 1.0936,
			"feet" => $this->earthRadius * 3.28084,
			"radians" => 1,
			"degrees" => $this->earthRadius / 111325,
		];
	}

	/**
	 * Boolean-contains returns True if the second geometry is completely contained by the first geometry.
	 * The interiors of both geometries must intersect and, the interior and boundary of the secondary (geometry b)
	 * must not intersect the exterior of the primary (geometry a).
	 * Boolean-contains returns the exact opposite result of the `@turf/boolean-within`.
	 *
	 * @name booleanContains
	 * @param {Geometry|Feature<any>} feature1 GeoJSON Feature or Geometry
	 * @param {Geometry|Feature<any>} feature2 GeoJSON Feature or Geometry
	 * @returns {boolean} true/false
	 * @example
	 * var line = turf.lineString([[1, 1], [1, 2], [1, 3], [1, 4]]);
	 * var point = turf.point([1, 2]);
	 *
	 * turf.booleanContains(line, point);
	 * //=true
	 */
	function booleanContains($feature1, $feature2) {
		$type1 = $this->getType($feature1);
		$type2 = $this->getType($feature2);
		$geom1 = $this->getGeom($feature1);
		$geom2 = $this->getGeom($feature2);
		$coords1 = $this->getCoords($feature1);
		$coords2 = $this->getCoords($feature2);

		switch ($type1) {
			case 'Point':
				switch ($type2) {
					case 'Point':
						return $this->compareCoords_2($coords1, $coords2);
					default:
						throw new InvalidArgumentException('feature2 ' . $type2 . ' geometry not supported');
				}
			case 'MultiPoint':
				switch ($type2) {
					case 'Point':
						return $this->isPointInMultiPoint_1($geom1, $geom2);
					case 'MultiPoint':
						return $this->isMultiPointInMultiPoint_1($geom1, $geom2);
					default:
						throw new InvalidArgumentException('feature2 ' . $type2 . ' geometry not supported');
				}
			case 'LineString':
				switch ($type2) {
					case 'Point':
						return $this->booleanPointOnLine($geom2, $geom1, (object) ['ignoreEndVertices' => true]);
			case 'LineString':
				return $this->isLineOnLine_2($geom1, $geom2);
			case 'MultiPoint':
				return $this->isMultiPointOnLine_1($geom1, $geom2);
			default:
				throw new InvalidArgumentException('feature2 ' . $type2 . ' geometry not supported');
		}
		case 'Polygon':
			switch ($type2) {
				case 'Point':
					return $this->booleanPointInPolygon($geom2, $geom1, (object) ['ignoreBoundary' => true ]);
				case 'LineString':
					return $this->isLineInPoly_2($geom1, $geom2);
				case 'Polygon':
					return $this->isPolyInPoly_2($geom1, $geom2);
				case 'MultiPoint':
					return $this->isMultiPointInPoly_1($geom1, $geom2);
				default:
					throw new InvalidArgumentException('feature2 ' . $type2 . ' geometry not supported');
			}
		default:
			throw new InvalidArgumentException('feature1 ' . $type1 . ' geometry not supported');
		}
	}

	function isMultiPointInPoly_1($polygon, $multiPoint) {
	    for ($i = 0; $i < count($multiPoint->coordinates); $i++) {
	        if (!$this->booleanPointInPolygon($multiPoint->coordinates[$i], $polygon, (object) ['ignoreBoundary' => true])) {
				return false;
			}
		}
		return true;
	}

	/**
	 * Is Polygon2 in Polygon1
	 * Only takes into account outer rings
	 *
	 * @private
	 * @param {Geometry|Feature<Polygon>} feature1 Polygon1
	 * @param {Geometry|Feature<Polygon>} feature2 Polygon2
	 * @returns {boolean} true/false
	 */
	function isPolyInPoly_2($feature1, $feature2) {
	    $poly1Bbox = $this->bbox($feature1);
	    $poly2Bbox = $this->bbox($feature2);

	    if (!$this->doBBoxOverlap_1($poly1Bbox, $poly2Bbox)) {
	        return false;
	    }
		for ($i = 0; $i < count($feature2->coordinates[0]); $i++) {
			if (!$this->booleanPointInPolygon($feature2->coordinates[0][$i], $feature1)) {
				return false;
			}
		}
	    return true;
	}

	function isLineInPoly_2($polygon, $linestring) {
	    $output = false;
	    $polyBbox = $this->bbox($polygon);
	    $lineBbox = $this->bbox($linestring);
	    if (!$this->doBBoxOverlap_1($polyBbox, $lineBbox)) {
	        return false;
	    }
		for ($i = 0; $i < count($linestring->coordinates) - 1; $i++) {
			$midPoint = $this->getMidpoint_1($linestring->coordinates[$i], $linestring->coordinates[$i + 1]);
	        if ($this->booleanPointInPolygon((object) ['type' => 'Point', 'coordinates' => $midPoint], $polygon, (object)['ignoreBoundary' => true ])) {
		        $output = true;
		        break;
	        }
	    }
	    return $output;
	}

	function getMidpoint_1($pair1, $pair2) {
	    return [($pair1[0] + $pair2[0]) / 2, ($pair1[1] + $pair2[1]) / 2];
	}

	function doBBoxOverlap_1($bbox1, $bbox2) {
	    if ($bbox1[0] > $bbox2[0]) return false;
	    if ($bbox1[2] < $bbox2[2]) return false;
	    if ($bbox1[1] > $bbox2[1]) return false;
	    if ($bbox1[3] < $bbox2[3]) return false;
	    return true;
	}

	// http://en.wikipedia.org/wiki/Even%E2%80%93odd_rule
	// modified from: https://github.com/substack/point-in-polygon/blob/master/index.js
	// which was modified from http://www.ecse.rpi.edu/Homepages/wrf/Research/Short_Notes/pnpoly.html

	/**
	 * Takes a {@link Point} and a {@link Polygon} or {@link MultiPolygon} and determines if the point resides inside the polygon. The polygon can
	 * be convex or concave. The function accounts for holes.
	 *
	 * @name booleanPointInPolygon
	 * @param {Coord} point input point
	 * @param {Feature<Polygon|MultiPolygon>} polygon input polygon or multipolygon
	 * @param {Object} [options={}] Optional parameters
	 * @param {boolean} [options.ignoreBoundary=false] True if polygon boundary should be ignored when determining if the point is inside the polygon otherwise false.
	 * @returns {boolean} `true` if the Point is inside the Polygon; `false` if the Point is not inside the Polygon
	 * @example
	 * var pt = turf.point([-77, 44]);
	 * var poly = turf.polygon([[
	 *   [-81, 41],
	 *   [-81, 47],
	 *   [-72, 47],
	 *   [-72, 41],
	 *   [-81, 41]
	 * ]]);
	 *
	 * turf.booleanPointInPolygon(pt, poly);
	 * //= true
	 */
	function booleanPointInPolygon($point, $polygon, $options = null) {
		// Optional parameters
		if (isset($options) && !is_object($options)) {
			throw new InvalidArgumentException('options is invalid');
		}
		$ignoreBoundary = isset($options) ? $options->ignoreBoundary : false;

	    // validation
	    if (!$point) throw new InvalidArgumentException('point is required');
	    if (!$polygon) throw new InvalidArgumentException('polygon is required');

	    $pt = $this->getCoord($point);
	    $polys = $this->getCoords($polygon);
	    $type = isset($polygon->geometry) ? $polygon->geometry->type : $polygon->type;
	    $bbox = isset($polygon->bbox) ? $polygon->bbox : null;

	    // Quick elimination if point is not inside bbox

	    if (isset($bbox) && $this->inBBox($pt, $bbox) === false) {
		    return false;
	    }

	    // normalize to multipolygon
	    if ($type === 'Polygon') $polys = [$polys];

	    for ($i = 0, $insidePoly = false; $i < count($polys) && !$insidePoly; $i++) {
			// check if it is in the outer ring first
			if ($this->inRing($pt, $polys[$i][0], $ignoreBoundary)) {
				$inHole = false;
				$k = 1;
				// check for the point in any of the holes
				while ($k < count($polys[$i]) && !$inHole) {
					if ($this->inRing($pt, $polys[$i][$k], !$ignoreBoundary)) {
						$inHole = true;
					}
					$k++;
				}
				if (!$inHole) $insidePoly = true;
			}
		}

	    return $insidePoly;
	}

	/**
	 * inRing
	 *
	 * @private
	 * @param {Array<number>} pt [x,y]
	 * @param {Array<Array<number>>} ring [[x,y], [x,y],..]
	 * @param {boolean} ignoreBoundary ignoreBoundary
	 * @returns {boolean} inRing
	 */
	function inRing($pt, $ring, $ignoreBoundary) {
		$isInside = false;

		if ($ring[0][0] === $ring[count($ring) - 1][0] && $ring[0][1] === $ring[count($ring) - 1][1]) $ring = array_splice($ring, 0, count($ring) - 1);

		for ($i = 0, $j = count($ring) - 1; $i < count($ring); $j = $i++) {
			$xi = $ring[$i][0]; $yi = $ring[$i][1];
	        $xj = $ring[$j][0]; $yj = $ring[$j][1];
	        $onBoundary = ($pt[1] * ($xi - $xj) + $yi * ($xj - $pt[0]) + $yj * ($pt[0] - $xi) === 0) &&
	                         (($xi - $pt[0]) * ($xj - $pt[0]) <= 0) && (($yi - $pt[1]) * ($yj - $pt[1]) <= 0);
	        if ($onBoundary) return !$ignoreBoundary;
	        $intersect = (($yi > $pt[1]) !== ($yj > $pt[1])) &&
	                        ($pt[0] < ($xj - $xi) * ($pt[1] - $yi) / ($yj - $yi) + $xi);
	        if ($intersect) $isInside = !$isInside;
	    }
	    return $isInside;
	}

	/**
	 * inBBox
	 *
	 * @private
	 * @param {Position} pt point [x,y]
	 * @param {BBox} bbox BBox [west, south, east, north]
	 * @returns {boolean} true/false if point is inside BBox
	 */
	function inBBox($pt, $bbox) {
		return $bbox[0] <= $pt[0] &&
		       $bbox[1] <= $pt[1] &&
		       $bbox[2] >= $pt[0] &&
		       $bbox[3] >= $pt[1];
	}

	function isMultiPointOnLine_1($lineString, $multiPoint) {
	    $haveFoundInteriorPoint = false;
	    for ($i = 0; $i < count($multiPoint->coordinates); $i++) {
	        if ($this->booleanPointOnLine($multiPoint->coordinates[$i], $lineString, (object) ['ignoreEndVertices' => true])) {
				$haveFoundInteriorPoint = true;
			}
			if (!$this->booleanPointOnLine($multiPoint->coordinates[$i], $lineString)) {
				return false;
			}
		}
		if ($haveFoundInteriorPoint) {
			return true;
		}
		return false;
	}

	function isLineOnLine_2($lineString1, $lineString2) {
	    $haveFoundInteriorPoint = false;
	    for ($i = 0; $i < count($lineString2->coordinates); $i++) {
	        if ($this->booleanPointOnLine((object)['type'=> 'Point', 'coordinates'=> $lineString2->coordinates[$i]], $lineString1, (object)['ignoreEndVertices' => true ])) {
				$haveFoundInteriorPoint = true;
			}
	        if (!$this->booleanPointOnLine((object)['type'=> 'Point', 'coordinates' => $lineString2->coordinates[$i]], $lineString1, (object)['ignoreEndVertices' => false ])) {
		        return false;
	        }
	    }
	    return $haveFoundInteriorPoint;
	}

	/**
	 * Returns true if a point is on a line. Accepts a optional parameter to ignore the start and end vertices of the linestring.
	 *
	 * @name booleanPointOnLine
	 * @param {Coord} pt GeoJSON Point
	 * @param {Feature<LineString>} line GeoJSON LineString
	 * @param {Object} [options={}] Optional parameters
	 * @param {boolean} [options.ignoreEndVertices=false] whether to ignore the start and end vertices.
	 * @returns {boolean} true/false
	 * @example
	 * var pt = turf.point([0, 0]);
	 * var line = turf.lineString([[-1, -1],[1, 1],[1.5, 2.2]]);
	 * var isPointOnLine = turf.booleanPointOnLine(pt, line);
	 * //=true
	 */
	function booleanPointOnLine($pt, $line, $options = null) {
		// Optional parameters
		$options = $options || (object) [];
	    $ignoreEndVertices = $options->ignoreEndVertices;
	    if (!is_object($options)) throw new InvalidArgumentException('invalid options');

	    // Validate input
	    if (!$pt) throw new InvalidArgumentException('pt is required');
	    if (!$line) throw new InvalidArgumentException('line is required');

	    // Normalize inputs
	    $ptCoords = $this->getCoord($pt);
	    $lineCoords = $this->getCoords($line);

	    // Main
	    for ($i = 0; $i < count($lineCoords) - 1; $i++) {
			$ignoreBoundary = false;
			if ($ignoreEndVertices) {
				if ($i === 0) $ignoreBoundary = 'start';
				if ($i === count($lineCoords) - 2) $ignoreBoundary = 'end';
				if ($i === 0 && $i + 1 === count($lineCoords) - 1) $ignoreBoundary = 'both';
			}
			if ($this->isPointOnLineSegment_1($lineCoords[$i], $lineCoords[$i + 1], $ptCoords, $ignoreBoundary)) return true;
	    }
	    return false;
	}

	// See http://stackoverflow.com/a/4833823/1979085
	/**
	 * @private
	 * @param {Position} lineSegmentStart coord pair of start of line
	 * @param {Position} lineSegmentEnd coord pair of end of line
	 * @param {Position} pt coord pair of point to check
	 * @param {boolean|string} excludeBoundary whether the point is allowed to fall on the line ends. If true which end to ignore.
	 * @returns {boolean} true/false
	 */
	function isPointOnLineSegment_1($lineSegmentStart, $lineSegmentEnd, $pt, $excludeBoundary) {
	    $x = $pt[0];
	    $y = $pt[1];
	    $x1 = $lineSegmentStart[0];
	    $y1 = $lineSegmentStart[1];
	    $x2 = $lineSegmentEnd[0];
	    $y2 = $lineSegmentEnd[1];
	    $dxc = $pt[0] - $x1;
	    $dyc = $pt[1] - $y1;
	    $dxl = $x2 - $x1;
	    $dyl = $y2 - $y1;
	    $cross = $dxc * $dyl - $dyc * $dxl;
	    if ($cross !== 0) {
	        return false;
	    }
		if (!$excludeBoundary) {
			if (abs($dxl) >= abs($dyl)) {
				return $dxl > 0 ? $x1 <= $x && $x <= $x2 : $x2 <= $x && $x <= $x1;
			}
			return $dyl > 0 ? $y1 <= $y && $y <= $y2 : $y2 <= $y && $y <= $y1;
		} else if ($excludeBoundary === 'start') {
			if (abs($dxl) >= abs($dyl)) {
				return $dxl > 0 ? $x1 < $x && $x <= $x2 : $x2 <= $x && $x < $x1;
			}
			return $dyl > 0 ? $y1 < $y && $y <= $y2 : $y2 <= $y && $y < $y1;
		} else if ($excludeBoundary === 'end') {
			if (abs($dxl) >= abs($dyl)) {
				return $dxl > 0 ? $x1 <= $x && $x < $x2 : $x2 < $x && $x <= $x1;
			}
			return $dyl > 0 ? $y1 <= $y && $y < $y2 : $y2 < $y && $y <= $y1;
		} else if ($excludeBoundary === 'both') {
			if (abs($dxl) >= abs($dyl)) {
				return $dxl > 0 ? $x1 < $x && $x < $x2 : $x2 < $x && $x < $x1;
			}
			return $dyl > 0 ? $y1 < $y && $y < $y2 : $y2 < $y && $y < $y1;
		}
	}

	/**
	 * compareCoords
	 *
	 * @private
	 * @param {Position} pair1 point [x,y]
	 * @param {Position} pair2 point [x,y]
	 * @returns {boolean} true/false if coord pairs match
	 */
	function compareCoords_1($pair1, $pair2) {
	    return $pair1[0] === $pair2[0] && $pair1[1] === $pair2[1];
	}

	function isPointInMultiPoint_1($multiPoint, $point) {
		$output = false;
		for ($i = 0; $i < count($multiPoint->coordinates); $i++) {
			if ($this->compareCoords_2($multiPoint->coordinates[$i], $point->coordinates)) {
				$output = true;
				break;
			}
		}
		return $output;
	}

	/**
	 * compareCoords
	 *
	 * @private
	 * @param {Position} pair1 point [x,y]
	 * @param {Position} pair2 point [x,y]
	 * @returns {boolean} true/false if coord pairs match
	 */
	function compareCoords_2($pair1, $pair2) {
	    return $pair1[0] === $pair2[0] && $pair1[1] === $pair2[1];
	}

	function isMultiPointInMultiPoint_1($multiPoint1, $multiPoint2) {
	    for ($i = 0; $i < count($multiPoint2->coordinates); $i++) {
	        $matchFound = false;
	        for ($i2 = 0; $i2 < count($multiPoint1->coordinates); $i2++) {
	            if ($this->compareCoords_2($multiPoint2->coordinates[$i], $multiPoint1->coordinates[$i2])) {
	                $matchFound = true;
	                break;
	            }
			}
			if (!$matchFound) {
				return false;
			}
		}
		return true;
	}

	/**
	 * Validate BBox
	 *
	 * @private
	 * @param {Array<number>} bbox BBox to validate
	 * @returns {void}
	 * @example
	 * validateBBox([-180, -40, 110, 50])
	 * //=OK
	 * validateBBox([-180, -40])
	 * //=Error
	 * validateBBox('Foo')
	 * //=Error
	 * validateBBox(5)
	 * //=Error
	 * validateBBox(null)
	 * //=Error
	 * validateBBox(undefined)
	 * //=Error
	 */
	function validateBBox($bbox) {
		if (!$bbox) throw new InvalidArgumentException('bbox is required');
		if (!is_array($bbox)) throw new InvalidArgumentException('bbox must be an Array');
		if (count($bbox) !== 4 && count($bbox) !== 6) throw new InvalidArgumentException('bbox must be an Array of 4 or 6 numbers');
		foreach ($bbox as $num) {
			if (!is_numeric($num)) throw new InvalidArgumentException('bbox must only contain numbers');
		}
	}

	/**
	 * Takes a GeoJSON Feature or FeatureCollection and truncates the precision of the geometry.
	 *
	 * @name truncate
	 * @param {GeoJSON} geojson any GeoJSON Feature, FeatureCollection, Geometry or GeometryCollection.
	 * @param {Object} [options={}] Optional parameters
	 * @param {number} [options.precision=6] coordinate decimal precision
	 * @param {number} [options.coordinates=3] maximum number of coordinates (primarly used to remove z coordinates)
	 * @param {boolean} [options.mutate=false] allows GeoJSON input to be mutated (significant performance increase if true)
	 * @returns {GeoJSON} layer with truncated geometry
	 * @example
	 * var point = turf.point([
	 *     70.46923055566859,
	 *     58.11088890802906,
	 *     1508
	 * ]);
	 * var options = {precision: 3, coordinates: 2};
	 * var truncated = turf.truncate(point, options);
	 * //=truncated.geometry.coordinates => [70.469, 58.111]
	 *
	 * //addToMap
	 * var addToMap = [truncated];
	 */
	function truncate($geojson, $options) {
		// Optional parameters
		$options = $options || (object) [];
	    if (!is_object($options)) throw new InvalidArgumentException('options is invalid');
	    $precision = $options->precision;
	    $coordinates = $options->coordinates;
	    $mutate = $options->mutate;

	    // default params
	    $precision = (!$precision || !is_numeric($precision)) ? 6 : $precision;
	    $coordinates = (!$coordinates || !is_numeric($coordinates)) ? 3 : $coordinates;

	    // validation
	    if (!$geojson) throw new InvalidArgumentException('<geojson> is required');
	    if (!is_numeric($precision)) throw new InvalidArgumentException('<precision> must be a number');
	    if (!is_numeric($coordinates)) throw new InvalidArgumentException('<coordinates> must be a number');

	    // prevent input mutation
	    //if (!$mutate) $geojson = JSON.parse(JSON.stringify(geojson));

	    $factor = pow(10, $precision);

	    // Truncate Coordinates
	    $this->coordEach($geojson, function (&$coords) use(&$factor, &$coordinates) {
		    $this->truncateCoords($coords, $factor, $coordinates);
	    });

	    return $geojson;
	}

	/**
	 * Truncate Coordinates - Mutates coordinates in place
	 *
	 * @private
	 * @param {Array<any>} coords Geometry Coordinates
	 * @param {number} factor rounding factor for coordinate decimal precision
	 * @param {number} coordinates maximum number of coordinates (primarly used to remove z coordinates)
	 * @returns {Array<any>} mutated coordinates
	 */
	function truncateCoords($coords, $factor, $coordinates) {
		// Remove extra coordinates (usually elevation coordinates and more)
		if (count($coords) > $coordinates) $coords = array_splice($coords, $coordinates, count($coords));

		// Truncate coordinate decimals
		for ($i = 0; $i < count($coords); $i++) {
			$coords[$i] = round($coords[$i] * $factor) / $factor;
		}
    return $coords;
}

	/**
	 * Validate Id
	 *
	 * @private
	 * @param {string|number} id Id to validate
	 * @returns {void}
	 * @example
	 * validateId([-180, -40, 110, 50])
	 * //=Error
	 * validateId([-180, -40])
	 * //=Error
	 * validateId('Foo')
	 * //=OK
	 * validateId(5)
	 * //=OK
	 * validateId(null)
	 * //=Error
	 * validateId(undefined)
	 * //=Error
	 */
	function validateId($id) {
		if (!$id) throw new InvalidArgumentException('id is required');
		if (!is_numeric($id) && !is_string($id)) {
			throw new InvalidArgumentException('id must be a number or a string');
		}
	}

	/**
	 * Get GeoJSON object's type, Geometry type is prioritize.
	 *
	 * @param {GeoJSON} geojson GeoJSON object
	 * @param {string} [name] name of the variable to display in error message
	 * @returns {string} GeoJSON type
	 * @example
	 * var point = {
	 *   "type": "Feature",
	 *   "properties": {},
	 *   "geometry": {
	 *     "type": "Point",
	 *     "coordinates": [110, 40]
	 *   }
	 * }
	 * var geom = turf.getType(point)
	 * //="Point"
	 */
	function getType($geojson, $name = "") {
		if (!$geojson) throw new InvalidArgumentException(($name || 'geojson') . ' is required');
		// GeoJSON Feature & GeometryCollection
		if ($geojson->geometry && $geojson->geometry->type) return $geojson->geometry->type;
		// GeoJSON Geometry & FeatureCollection
		if ($geojson->type) return $geojson->type;
		throw new InvalidArgumentException(($name || 'geojson') . ' is invalid');
	}

	/**
	 * Converts an angle in degrees to radians
	 *
	 * @name degreesToRadians
	 * @param {number} degrees angle between 0 and 360 degrees
	 * @returns {number} angle in radians
	 */
	function degreesToRadians($degrees) {
		if (!isset($degrees)) throw new InvalidArgumentException('degrees is required');

		$radians = fmod($degrees, 360);
		return $radians * M_PI / 180;
	}

	/**
	 * Convert a distance measurement (assuming a spherical Earth) from radians to a more friendly unit.
	 * Valid units: miles, nauticalmiles, inches, yards, meters, metres, kilometers, centimeters, feet
	 *
	 * @name radiansToLength
	 * @param {number} radians in radians across the sphere
	 * @param {string} [units='kilometers'] can be degrees, radians, miles, or kilometers inches, yards, metres, meters, kilometres, kilometers.
	 * @returns {number} distance
	 */
	function radiansToLength($radians, $units) {
		if (!isset($radians)) throw new InvalidArgumentException('radians is required');

		if ($units && ! is_string($units)) throw new InvalidArgumentException('units must be a string');
	    $factor = $this->factors[isset($units) ? $units : 'kilometers'];
	    if (!$factor) throw new InvalidArgumentException($units . ' units is invalid');
	    return $radians * $factor;
	}

	//http://en.wikipedia.org/wiki/Haversine_formula
	//http://www.movable-type.co.uk/scripts/latlong.html

	/**
	 * Calculates the distance between two {@link Point|points} in degrees, radians,
	 * miles, or kilometers. This uses the
	 * [Haversine formula](http://en.wikipedia.org/wiki/Haversine_formula)
	 * to account for global curvature.
	 *
	 * @name distance
	 * @param {Coord} from origin point
	 * @param {Coord} to destination point
	 * @param {Object} [options={}] Optional parameters
	 * @param {string} [options.units='kilometers'] can be degrees, radians, miles, or kilometers
	 * @returns {number} distance between the two points
	 * @example
	 * var from = turf.point([-75.343, 39.984]);
	 * var to = turf.point([-75.534, 39.123]);
	 * var options = {units: 'miles'};
	 *
	 * var distance = turf.distance(from, to, options);
	 *
	 * //addToMap
	 * var addToMap = [from, to];
	 * from.properties.distance = distance;
	 * to.properties.distance = distance;
	 */
	function distance($from, $to, $options) {
		// Optional parameters
		if (isset($options) && !is_object($options)) {
			throw new InvalidArgumentException('options is invalid');
		}
		$units = isset($options->units) ? $options->units : '';

		$coordinates1 = $this->getCoord($from);
		$coordinates2 = $this->getCoord($to);
		$dLat = $this->degreesToRadians(($coordinates2[1] - $coordinates1[1]));
		$dLon = $this->degreesToRadians(($coordinates2[0] - $coordinates1[0]));
		$lat1 = $this->degreesToRadians($coordinates1[1]);
		$lat2 = $this->degreesToRadians($coordinates2[1]);

		$a = pow(sin($dLat / 2), 2) +
			pow(sin($dLon / 2), 2) * cos($lat1) * cos($lat2);

		return $this->radiansToLength(2 * atan2(sqrt($a), sqrt(1 - $a)), $units);
	}

	/**
	 * Checks if coordinates contains a number
	 *
	 * @name containsNumber
	 * @param {Array<any>} coordinates GeoJSON Coordinates
	 * @returns {boolean} true if Array contains a number
	 */
	function containsNumber($coordinates) {
		if (count($coordinates) > 1 && is_numeric($coordinates[0]) && is_numeric($coordinates[1])) {
			return true;
		}

		if (is_array($coordinates[0]) && count($coordinates[0])) {
			return $this->containsNumber($coordinates[0]);
		}
		throw new UnexpectedValueException('coordinates must only contain numbers');
	}

	/**
	 * Get Geometry from Feature or Geometry Object
	 *
	 * @param {Feature|Geometry} geojson GeoJSON Feature or Geometry Object
	 * @returns {Geometry|null} GeoJSON Geometry Object
	 * @throws {Error} if geojson is not a Feature or Geometry Object
	 * @example
	 * var point = {
	 *   "type": "Feature",
	 *   "properties": {},
	 *   "geometry": {
	 *     "type": "Point",
	 *     "coordinates": [110, 40]
	 *   }
	 * }
	 * var geom = turf.getGeom(point)
	 * //={"type": "Point", "coordinates": [110, 40]}
	 */
	function getGeom($geojson) {
		if (!$geojson) throw new InvalidArgumentException('geojson is required');
		if ($geojson->geometry) return $geojson->geometry;
		if ($geojson->coordinates || $geojson->geometries) return $geojson;
		throw new UnexpectedValueException('geojson must be a valid Feature or Geometry Object');
	}

	/**
	 * Unwrap a coordinate from a Point Feature, Geometry or a single coordinate.
	 *
	 * @name getCoord
	 * @param {Array<number>|Geometry<Point>|Feature<Point>} obj Object
	 * @returns {Array<number>} coordinates
	 * @example
	 * var pt = turf.point([10, 10]);
	 *
	 * var coord = turf.getCoord(pt);
	 * //= [10, 10]
	 */
	function getCoord($obj) {
		if (!$obj) throw new InvalidArgumentException('obj is required');

		$coordinates = $this->getCoords($obj);

		// getCoord() must contain at least two numbers (Point)
		if (count($coordinates) > 1 && is_numeric($coordinates[0]) && is_numeric($coordinates[1])) {
			return $coordinates;
		} else {
			throw new UnexpectedValueException('Coordinate is not a valid Point');
		}
	}

	/**
	 * Unwrap coordinates from a Feature, Geometry Object or an Array of numbers
	 *
	 * @name getCoords
	 * @param {Array<number>|Geometry|Feature} obj Object
	 * @returns {Array<number>} coordinates
	 * @example
	 * var poly = turf.polygon([[[119.32, -8.7], [119.55, -8.69], [119.51, -8.54], [119.32, -8.7]]]);
	 *
	 * var coord = turf.getCoords(poly);
	 * //= [[[119.32, -8.7], [119.55, -8.69], [119.51, -8.54], [119.32, -8.7]]]
	 */
	function getCoords($obj) {
		if (!$obj) throw new InvalidArgumentException('obj is required');
		$coordinates = '';

		// Array of numbers
		if (is_array($obj)) {
			$coordinates = $obj;

			// Geometry Object
		} else if (isset($obj->coordinates)) {
			$coordinates = $obj->coordinates;

			// Feature
		} else if (isset($obj->geometry) && isset($obj->geometry->coordinates)) {
			$coordinates = $obj->geometry->coordinates;
		}
		// Checks if coordinates contains a number
		if ($coordinates) {
			$this->containsNumber($coordinates);
			return $coordinates;
		}
		throw new UnexpectedValueException('No valid coordinates');
	}

	/**
	 * Removes redundant coordinates from any GeoJSON Geometry.
	 *
	 * @name cleanCoords
	 * @param {Geometry|Feature} geojson Feature or Geometry
	 * @param {Object} [options={}] Optional parameters
	 * @param {boolean} [options.mutate=false] allows GeoJSON input to be mutated
	 * @returns {Geometry|Feature} the cleaned input Feature/Geometry
	 * @example
	 * var line = turf.lineString([[0, 0], [0, 2], [0, 5], [0, 8], [0, 8], [0, 10]]);
	 * var multiPoint = turf.multiPoint([[0, 0], [0, 0], [2, 2]]);
	 *
	 * turf.cleanCoords(line).geometry.coordinates;
	 * //= [[0, 0], [0, 10]]
	 *
	 * turf.cleanCoords(multiPoint).geometry.coordinates;
	 * //= [[0, 0], [2, 2]]
	 */
	function cleanCoords($geojson, $options) {
		// Backwards compatible with v4.0
		$mutate = is_object($options) ? $options->mutate : $options;
	    if (!$geojson) throw new InvalidArgumentException('geojson is required');
	    $type = $this->getType($geojson);

	    // Store new "clean" points in this Array
	    $newCoords = [];

	    switch ($type) {
		    case 'LineString':
			    $newCoords = $this->cleanLine($geojson);
			    break;
		    case 'MultiLineString':
		    case 'Polygon':
			    $lines = $this->getCoords($geojson);
			    foreach ($lines as $line) {
				    $newCoords[$this->cleanLine($line)];
			    };
	            break;
		    case 'MultiPolygon':
			    $polygons = $this->getCoords($geojson);
			    foreach ($polygons as $polygons__1) {
				    $polyPoints = [];
				    foreach ($polygons__1 as $ring) {
					    $polyPoints[$this->cleanLine($ring)];
				    }
	                $newCoords[$polyPoints];
			    }
	            break;
		    case 'Point':
			    return $geojson;
		    case 'MultiPoint':
			    $existing = (object) [];
			    $coords = $this->getCoords($geojson);
			    foreach ($coords as $coord) {
				    $key = join('-', $coord);
				    if (!property_exists($existing, $key)) {
					    $newCoords[$coord];
					    $existing[$key] = true;
				    }
			    }
	            break;
		    default:
			    throw new UnexpectedValueException($type . ' geometry not supported');
	    }

	    // Support input mutation
	    if ($geojson->coordinates) {
		    if ($mutate === true) {
			    $geojson->coordinates = $newCoords;
			    return $geojson;
		    }
		    return (object) ['type' => $type, 'coordinates' => 'newCoords'];
	    } else {
		    if ($mutate === true) {
			    $geojson->geometry->coordinates = $newCoords;
			    return $geojson;
		    }
		    //return $this->feature((object) ['type' => $type, 'coordinates' => $newCoords], $geojson->properties, $geojson->bbox, $geojson->id);
		    return $this->feature((object) ['type' => $type, 'coordinates' => $newCoords], $geojson->properties, $geojson->bbox);
	    }
	}

	/**
	 * Compares two points and returns if they are equals
	 *
	 * @private
	 * @param {Position} pt1 point
	 * @param {Position} pt2 point
	 * @returns {boolean} true if they are equals
	 */
	function equals($pt1, $pt2) {
		return $pt1[0] === $pt2[0] && $pt1[1] === $pt2[1];
	}

	/**
	 * Returns if `point` is on the segment between `start` and `end`.
	 * Borrowed from `@turf/boolean-point-on-line` to speed up the evaluation (instead of using the module as dependency)
	 *
	 * @private
	 * @param {Position} start coord pair of start of line
	 * @param {Position} end coord pair of end of line
	 * @param {Position} point coord pair of point to check
	 * @returns {boolean} true/false
	 */
	function isPointOnLineSegment($start, $end, $point__1) {
	    $x = $point__1[0]; $y = $point__1[1];
	    $startX = $start[0]; $startY = $start[1];
	    $endX = $end[0]; $endY = $end[1];

	    $dxc = $x - $startX;
	    $dyc = $y - $startY;
	    $dxl = $endX - $startX;
	    $dyl = $endY - $startY;
	    $cross = $dxc * $dyl - $dyc * $dxl;

	    if ($cross !== 0) return false;

	    else if (abs($dxl) >= abs($dyl)) return $dxl > 0 ? $startX <= $x && $x <= $endX : $endX <= $x && $x <= $startX;
	    else return $dyl > 0 ? $startY <= $y && $y <= $endY : $endY <= $y && $y <= $startY;
	}

	/**
	 * Clean Coords
	 *
	 * @private
	 * @param {Array<number>|LineString} line Line
	 * @returns {Array<number>} Cleaned coordinates
	 */
	function cleanLine($line) {
		$points__1 = $this->getCoords($line);
	    // handle "clean" segment
	    if (count($points__1) === 2 && !$this->equals($points__1[0], $points__1[1])) return $points__1;

	    $newPoints = [];
		$nextPoint = "";
	    $secondToLast = count($points__1) - 1;

	    $newPoints[$points__1[0]];
	    for ($i = 1; $i < $secondToLast; $i++) {
			$prevPoint = $points__1[$i - 1];
	        $point__1 = $points__1[$i];
	        $nextPoint = $points__1[$i + 1];

	        if (!$this->isPointOnLineSegment($prevPoint, $nextPoint, $point__1)) {
					$newPoints[$point__1];
	        }
	    }
	    $newPoints[$nextPoint];
	    return $newPoints;
	}


	/**
	 * Creates a {@link Polygon} {@link Feature} from an Array of LinearRings.
	 *
	 * @name polygon
	 * @param {Array<Array<Array<number>>>} coordinates an array of LinearRings
	 * @param {Object} [properties={}] an Object of key-value pairs to add as properties
	 * @param {Object} [options={}] Optional Parameters
	 * @param {Array<number>} [options.bbox] Bounding Box Array [west, south, east, north] associated with the Feature
	 * @param {string|number} [options.id] Identifier associated with the Feature
	 * @returns {Feature<Polygon>} Polygon Feature
	 * @example
	 * var polygon = turf.polygon([[[-5, 52], [-4, 56], [-2, 51], [-7, 54], [-5, 52]]], { name: 'poly1' });
	 *
	 * //=polygon
	 */
	function polygon($coordinates, $properties = null, $options = null) {
		if (!$coordinates) throw new InvalidArgumentException('coordinates is required');

		for ($i = 0; $i < count($coordinates); $i++) {
			$ring = $coordinates[$i];
			$ringLength = count($ring);
			if ($ringLength < 4) {
				throw new InvalidArgumentException('Each LinearRing of a Polygon must have 4 or more Positions.');
			}
			for ($j = 0; $j < count($ring[$ringLength - 1]); $j++) {
				// Check if first point of Polygon contains two numbers
				if ($i === 0 && $j === 0 && !is_numeric($ring[0][0]) || !is_numeric($ring[0][1])) throw new InvalidArgumentException('coordinates must contain numbers');
				if ($ring[$ringLength - 1][$j] !== $ring[0][$j]) {
					throw new InvalidArgumentException('First and last Position are not equivalent.');
				}
			}
		}

		$feature = (object) [
			"type" => 'Polygon',
			"coordinates" => $coordinates
		];

		return $this->feature($feature, $properties, $options);
	}

	/**
	 * Wraps a GeoJSON {@link Geometry} in a GeoJSON {@link Feature}.
	 *
	 * @name feature
	 * @param {Geometry} geometry input geometry
	 * @param {Object} [properties={}] an Object of key-value pairs to add as properties
	 * @param {Object} [options={}] Optional Parameters
	 * @param {Array<number>} [options.bbox] Bounding Box Array [west, south, east, north] associated with the Feature
	 * @param {string|number} [options.id] Identifier associated with the Feature
	 * @returns {Feature} a GeoJSON Feature
	 * @example
	 * var geometry = {
	 *   "type": "Point",
	 *   "coordinates": [110, 50]
	 * };
	 *
	 * var feature = turf.feature(geometry);
	 *
	 * //=feature
	 */
	function feature($geometry, $properties, $options = null) {
		// Optional Parameters
		$bbox = '';
		$id = '';
		if (isset($options)) {
			if (!is_object($options)) {
				throw new InvalidArgumentException('options is invalid');
			}
			$bbox = $options->bbox;
			$id = $options->id;
		}

		// Validation
		if (!$geometry) throw new InvalidArgumentException('geometry is required');
		if ($properties && ! is_object($properties->constructor)) throw new InvalidArgumentException('properties must be an Object');
		if ($bbox) $this->validateBBox($bbox);
		if ($id) $this->validateId($id);

		// Main
		$feat = (object) ['type' => 'Feature'];
		if ($bbox) $feat->id = $id;
		if ($id) $feat->bbox = $bbox;
		$feat->properties = $properties || (object) [];
		$feat->geometry = $geometry;

		return $feat;
	}

	/**
	 * Takes a set of features, calculates the bbox of all input features, and returns a bounding box.
	 *
	 * @name bbox
	 * @param {GeoJSON} geojson any GeoJSON object
	 * @returns {BBox} bbox extent in [minX, minY, maxX, maxY] order
	 * @example
	 * var line = turf.lineString([[-74, 40], [-78, 42], [-82, 35]]);
	 * var bbox = turf.bbox(line);
	 * var bboxPolygon = turf.bboxPolygon(bbox);
	 *
	 * //addToMap
	 * var addToMap = [line, bboxPolygon]
	 */
	function bbox($geojson) {
		$BBox = [INF, INF, -INF, -INF];

		$this->coordEach($geojson, function ($coord) use (&$BBox) {
			if ($BBox[0] > $coord[0]) $BBox[0] = $coord[0];
			if ($BBox[1] > $coord[1]) $BBox[1] = $coord[1];
			if ($BBox[2] < $coord[0]) $BBox[2] = $coord[0];
			if ($BBox[3] < $coord[1]) $BBox[3] = $coord[1];
		});

		return $BBox;
	}

	/**
	 * Callback for coordEach
	 *
	 * @callback coordEachCallback
	 * @param {Array<number>} currentCoord The current coordinate being processed.
	 * @param {number} coordIndex The current index of the coordinate being processed.
	 * @param {number} featureIndex The current index of the Feature being processed.
	 * @param {number} multiFeatureIndex The current index of the Multi-Feature being processed.
	 * @param {number} geometryIndex The current index of the Geometry being processed.
	 */

	/**
	 * Iterate over coordinates in any GeoJSON object, similar to Array.forEach()
	 *
	 * @name coordEach
	 * @param {FeatureCollection|Feature|Geometry} geojson any GeoJSON object
	 * @param {Function} callback a method that takes (currentCoord, coordIndex, featureIndex, multiFeatureIndex)
	 * @param {boolean} [excludeWrapCoord=false] whether or not to include the final coordinate of LinearRings that wraps the ring in its iteration.
	 * @example
	 * var features = turf.featureCollection([
	 *   turf.point([26, 37], {"foo": "bar"}),
	 *   turf.point([36, 53], {"hello": "world"})
	 * ]);
	 *
	 * turf.coordEach(features, function (currentCoord, coordIndex, featureIndex, multiFeatureIndex, geometryIndex) {
	 *   //=currentCoord
	 *   //=coordIndex
	 *   //=featureIndex
	 *   //=multiFeatureIndex
	 *   //=geometryIndex
	 * });
	 */
	function coordEach($geojson, $callback, $excludeWrapCoord = null) {
		// Handles null Geometry -- Skips this GeoJSON
		if ($geojson === null) return;
		$coordIndex = 0;
		$type = $geojson->type;
		$isFeatureCollection = $type === 'FeatureCollection';
		$isFeature = $type === 'Feature';
		$stop = $isFeatureCollection ? count($geojson->features) : 1;

		// This logic may look a little weird. The reason why it is that way
		// is because it's trying to be fast. GeoJSON supports multiple kinds
		// of objects at its root: FeatureCollection, Features, Geometries.
		// This function has the responsibility of handling all of them, and that
		// means that some of the `for` loops you see below actually just don't apply
		// to certain inputs. For instance, if you give this just a
		// Point geometry, then both loops are short-circuited and all we do
		// is gradually rename the input until it's called 'geometry'.
		//
		// This also aims to allocate as few resources as possible: just a
		// few numbers and booleans, rather than any temporary arrays as would
		// be required with the normalization approach.
		for ($featureIndex = 0; $featureIndex < $stop; $featureIndex++) {
			$geometryMaybeCollection = ($isFeatureCollection ? $geojson->features[$featureIndex]->geometry :
				($isFeature ? $geojson->geometry : $geojson));
			$isGeometryCollection = ($geometryMaybeCollection) ? $geometryMaybeCollection->type === 'GeometryCollection' : false;
			$stopG = $isGeometryCollection ? count($geometryMaybeCollection->geometries) : 1;

			for ($geomIndex = 0; $geomIndex < $stopG; $geomIndex++) {
				$multiFeatureIndex = 0;
				$geometryIndex = 0;
				$geometry__1 = $isGeometryCollection ?
					$geometryMaybeCollection->geometries[$geomIndex] : $geometryMaybeCollection;

				// Handles null Geometry -- Skips this geometry
				if ($geometry__1 === null) continue;
				$coords = $geometry__1->coordinates;
				$geomType = $geometry__1->type;

				$wrapShrink = ($excludeWrapCoord && ($geomType === 'Polygon' || $geomType === 'MultiPolygon')) ? 1 : 0;

				switch ($geomType) {
					case null:
						break;
					case 'Point':
						$callback($coords, $coordIndex, $featureIndex, $multiFeatureIndex, $geometryIndex);
						$coordIndex++;
						$multiFeatureIndex++;
						break;
					case 'LineString':
					case 'MultiPoint':
						for ($j = 0; $j < count($coords); $j++) {
							$callback($coords[$j], $coordIndex, $featureIndex, $multiFeatureIndex, $geometryIndex);
							$coordIndex++;
							if ($geomType === 'MultiPoint') $multiFeatureIndex++;
						}
						if ($geomType === 'LineString') $multiFeatureIndex++;
						break;
					case 'Polygon':
					case 'MultiLineString':
						for ($j = 0; $j < count($coords); $j++) {
							for ($k = 0; $k < count($coords[$j]) - $wrapShrink; $k++) {
								$callback($coords[$j][$k], $coordIndex, $featureIndex, $multiFeatureIndex, $geometryIndex);
								$coordIndex++;
							}
							if ($geomType === 'MultiLineString') $multiFeatureIndex++;
							if ($geomType === 'Polygon') $geometryIndex++;
						}
						if ($geomType === 'Polygon') $multiFeatureIndex++;
						break;
					case 'MultiPolygon':
						for ($j = 0; $j < count($coords); $j++) {
							if ($geomType === 'MultiPolygon') $geometryIndex = 0;
							for ($k = 0; $k < count($coords[$j]); $k++) {
								for ($l = 0; $l < count($coords[$j][$k]) - $wrapShrink; $l++) {
									$callback($coords[$j][$k][$l], $coordIndex, $featureIndex, $multiFeatureIndex, $geometryIndex);
									$coordIndex++;
								}
								$geometryIndex++;
							}
							$multiFeatureIndex++;
						}
						break;
					case 'GeometryCollection':
						for ($j = 0; $j < count($geometry__1->geometries); $j++)
						$this->coordEach($geometry__1->geometries[$j], $callback, $excludeWrapCoord);
						break;
					default:
						throw new InvalidArgumentException('Unknown Geometry Type');
				}
			}
		}
	}

	// depend on jsts for now http://bjornharrtell.github.io/jsts/
	/**
	 * Takes two {@link Polygon|polygons} and finds their intersection. If they share a border, returns the border; if they don't intersect, returns undefined.
	 *
	 * @name intersect
	 * @param {Feature<Polygon>} poly1 the first polygon
	 * @param {Feature<Polygon>} poly2 the second polygon
	 * @returns {Feature|null} returns a feature representing the point(s) they share (in case of a {@link Point}  or {@link MultiPoint}), the borders they share (in case of a {@link LineString} or a {@link MultiLineString}), the area they share (in case of {@link Polygon} or {@link MultiPolygon}). If they do not share any point, returns `null`.
	 * @example
	 * var poly1 = turf.polygon([[
	 *   [-122.801742, 45.48565],
	 *   [-122.801742, 45.60491],
	 *   [-122.584762, 45.60491],
	 *   [-122.584762, 45.48565],
	 *   [-122.801742, 45.48565]
	 * ]]);
	 *
	 * var poly2 = turf.polygon([[
	 *   [-122.520217, 45.535693],
	 *   [-122.64038, 45.553967],
	 *   [-122.720031, 45.526554],
	 *   [-122.669906, 45.507309],
	 *   [-122.723464, 45.446643],
	 *   [-122.532577, 45.408574],
	 *   [-122.487258, 45.477466],
	 *   [-122.520217, 45.535693]
	 * ]]);
	 *
	 * var intersection = turf.intersect(poly1, poly2);
	 *
	 * //addToMap
	 * var addToMap = [poly1, poly2, intersection];
	 */
	function intersect_2($poly1, $poly2) {
	    $geom1 = $this->getGeom($poly1);
	    $geom2 = $this->getGeom($poly2);

	    // Return null if geometry is too narrow in coordinate precision
	    // fixes topology errors with JSTS
	    // https://github.com/Turfjs/turf/issues/463
	    // https://github.com/Turfjs/turf/pull/1004
	    if (count($this->cleanCoords($this->truncate($geom2, (object) ['precision' => 4]))->coordinates[0]) < 4) return null;
	    if (count($this->cleanCoords($this->truncate($geom1, (object) ['precision' => 4]))->coordinates[0]) < 4) return null;

	    /*$reader = new GeoJSONReader();
	    var a = reader.read(truncate(geom1));
	    var b = reader.read(truncate(geom2));
	    var intersection = OverlayOp.intersection(a, b);

	    // https://github.com/Turfjs/turf/issues/951
	    if (intersection.isEmpty()) return null;

	    var writer = new GeoJSONWriter();
	    var geom = writer.write(intersection);
	    return feature(geom);*/
	}

	/**
	 * Creates a square grid from a bounding box, {@link Feature} or {@link FeatureCollection}.
	 *
	 * @name squareGrid
	 * @param {Array<number>} bbox extent in [minX, minY, maxX, maxY] order
	 * @param {number} cellSide of each cell, in units
	 * @param {Object} [options={}] Optional parameters
	 * @param {string} [options.units='kilometers'] used in calculating cellSide, can be degrees, radians, miles, or kilometers
	 * @param {Feature<Polygon|MultiPolygon>} [options.mask] if passed a Polygon or MultiPolygon, the grid Points will be created only inside it
	 * @param {Object} [options.properties={}] passed to each point of the grid
	 * @returns {FeatureCollection<Polygon>} grid a grid of polygons
	 * @example
	 * var bbox = [-95, 30 ,-85, 40];
	 * var cellSide = 50;
	 * var options = {units: 'miles'};
	 *
	 * var squareGrid = turf.squareGrid(bbox, cellSide, options);
	 *
	 * //addToMap
	 * var addToMap = [squareGrid]
	 */
	function squareGrid($bbox, $cellSide, $options) {
		// Optional parameters
		if (isset($options) && !is_object($options)) {
			throw new InvalidArgumentException('options is invalid');
		}
		$properties = isset($options->properties) ? $options->properties : '';
		$mask = isset($options->mask) ? $options->mask : '';

		// Containers
		$results = [];

		// Input Validation
		if ($cellSide === null) throw new InvalidArgumentException('cellSide is required');
		if (!is_numeric($cellSide)) throw new InvalidArgumentException('cellSide is invalid');
		if (!$bbox) throw new InvalidArgumentException('bbox is required');
		if (!is_array($bbox)) throw new InvalidArgumentException('bbox must be array');
		if (count($bbox) !== 4) throw new InvalidArgumentException('bbox must contain 4 numbers');
		if ($mask && !in_array($this->getType($mask), ['Polygon', 'MultiPolygon'])) throw new InvalidArgumentException('options.mask must be a (Multi)Polygon');

		$west = $bbox[0];
		$south = $bbox[1];
		$east = $bbox[2];
		$north = $bbox[3];

		$xFraction = $cellSide / ($this->distance([$west, $south], [$east, $south], $options));
		$cellWidth = $xFraction * ($east - $west);
		$yFraction = $cellSide / ($this->distance([$west, $south], [$west, $north], $options));
		$cellHeight = $yFraction * ($north - $south);

		// rows & columns
		$bboxWidth = ($east - $west);
		$bboxHeight = ($north - $south);
		$columns = floor($bboxWidth / $cellWidth);
		$rows = floor($bboxHeight / $cellHeight);

		// adjust origin of the grid
		$deltaX = ($bboxWidth - $columns * $cellWidth) / 2;
		$deltaY = ($bboxHeight - $rows * $cellHeight) / 2;

		// iterate over columns & rows
		$currentX = $west + $deltaX;
		for ($column = 0; $column < $columns; $column++) {
			$currentY = $south + $deltaY;
			for ($row = 0; $row < $rows; $row++) {
				$cellPoly = $this->polygon([[
				[$currentX, $currentY],
				[$currentX, $currentY + $cellHeight],
				[$currentX + $cellWidth, $currentY + $cellHeight],
				[$currentX + $cellWidth, $currentY],
				[$currentX, $currentY]
				]], $properties);
				if ($mask) {
					// This involves a lot of convert work from Turf JS
					// since we do not use mask, ignore it for now.
					//if ($this->intersect_2($mask, $cellPoly)) $results[$cellPoly];
				} else {
					$results[] = $cellPoly;
				}

				$currentY += $cellHeight;
			}

			$currentX += $cellWidth;
		}

		return $this->featureCollection($results);
	}

	function squareGridInfoOnly($bbox, $cellSide, $options) {
		// Optional parameters
		if (isset($options) && !is_object($options)) {
			throw new InvalidArgumentException('options is invalid');
		}
		$properties = isset($options->properties) ? $options->properties : '';
		$mask = isset($options->mask) ? $options->mask : '';

		// Containers
		$results = [];

		// Input Validation
		if ($cellSide === null) throw new InvalidArgumentException('cellSide is required');
		if (!is_numeric($cellSide)) throw new InvalidArgumentException('cellSide is invalid');
		if (!$bbox) throw new InvalidArgumentException('bbox is required');
		if (!is_array($bbox)) throw new InvalidArgumentException('bbox must be array');
		if (count($bbox) !== 4) throw new InvalidArgumentException('bbox must contain 4 numbers');
		if ($mask && !in_array($this->getType($mask), ['Polygon', 'MultiPolygon'])) throw new InvalidArgumentException('options.mask must be a (Multi)Polygon');

		$west = $bbox[0];
		$south = $bbox[1];
		$east = $bbox[2];
		$north = $bbox[3];

		$xFraction = $cellSide / ($this->distance([$west, $south], [$east, $south], $options));
		$cellWidth = $xFraction * ($east - $west);
		$yFraction = $cellSide / ($this->distance([$west, $south], [$west, $north], $options));
		$cellHeight = $yFraction * ($north - $south);

		// rows & columns
		$bboxWidth = ($east - $west);
		$bboxHeight = ($north - $south);
		$columns = floor($bboxWidth / $cellWidth);
		$rows = floor($bboxHeight / $cellHeight);

		// adjust origin of the grid
		$deltaX = ($bboxWidth - $columns * $cellWidth) / 2;
		$deltaY = ($bboxHeight - $rows * $cellHeight) / 2;

		return [
			'columns' => $columns,
			'rows' => $rows,
			'deltaX' => $deltaX,
			'deltaY' => $deltaY,
			'west' => $west,
			'south' => $south,
			'cellHeight' => $cellHeight,
			'cellWidth' => $cellWidth
		];
	}

	function squareGridByChunks($from, $chunks, $squareGridInfo) {
		$west = $squareGridInfo['west'];
		$south = $squareGridInfo['south'];
		$deltaX = $squareGridInfo['deltaX'];
		$deltaY = $squareGridInfo['deltaY'];
		$cellHeight = $squareGridInfo['cellHeight'];
		$cellWidth = $squareGridInfo['cellWidth'];
		$columns = $squareGridInfo['columns'];
		$rows = $squareGridInfo['rows'];
		$results = [];
		$counter = 0;

		// iterate over columns & rows
		$fromColumn = floor( $from / $rows );
		$currentX = ($west + $deltaX) + ($fromColumn * $cellWidth);
		for ($column = $fromColumn; $column < $columns; $column++) {
			$currentY = ($south + $deltaY) + ($cellHeight * ($from % $rows));
			for ($row = $from % $rows; $row <= $rows; $row++) {
				$cellPoly = $this->polygon([[
					[$currentX, $currentY],
					[$currentX, $currentY + $cellHeight],
					[$currentX + $cellWidth, $currentY + $cellHeight],
					[$currentX + $cellWidth, $currentY],
					[$currentX, $currentY]
				]]);

				if ( $counter == $chunks || ($from++ == $rows * $columns)) {
					break 2;
				}

				$results[] = $cellPoly;
				$currentY += $cellHeight;
				$counter++;
			}
			$currentX += $cellWidth;
		}

		return $results;
	}

	/**
	 * Takes one or more {@link Feature|Features} and creates a {@link FeatureCollection}.
	 *
	 * @name featureCollection
	 * @param {Feature[]} features input features
	 * @param {Object} [options={}] Optional Parameters
	 * @param {Array<number>} [options.bbox] Bounding Box Array [west, south, east, north] associated with the Feature
	 * @param {string|number} [options.id] Identifier associated with the Feature
	 * @returns {FeatureCollection} FeatureCollection of Features
	 * @example
	 * var locationA = turf.point([-75.343, 39.984], {name: 'Location A'});
	 * var locationB = turf.point([-75.833, 39.284], {name: 'Location B'});
	 * var locationC = turf.point([-75.534, 39.123], {name: 'Location C'});
	 *
	 * var collection = turf.featureCollection([
	 *   locationA,
	 *   locationB,
	 *   locationC
	 * ]);
	 *
	 * //=collection
	 */
	function featureCollection($features, $options = null) {
		// Optional Parameters
		$options = $options || (object)[];
	    if (!is_object($options)) throw new InvalidArgumentException('options is invalid');
	    $bbox = $options->bbox;
	    $id = $options->id;

	    // Validation
	    if (!$features) throw new InvalidArgumentException('No features passed');
	    if (!is_array($features)) throw new InvalidArgumentException('features must be an Array');
	    if ($bbox) $this->validateBBox($bbox);
	    if ($id) $this->validateId($id);

	    // Main
	    $fc = (object) ['type' => 'FeatureCollection'];
	    if ($id) $fc->id = $id;
	    if ($bbox) $fc->bbox = $bbox;
	    $fc->features = $features;

	    return $fc;
	}
}