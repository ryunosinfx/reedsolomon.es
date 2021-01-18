/*
 * Copyright 2007 ZXing authors
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

// package com.google.zxing.common.reedsolomon;

/**
 * <p>This class contains utility methods for performing mathematical operations over
 * the Galois Fields. Operations use a given primitive polynomial in calculations.</p>
 *
 * <p>Throughout this package, elements of the GF are represented as an {@code int}
 * for convenience and speed (but at the cost of memory).
 * </p>
 *
 * @author Sean Owen
 * @author David Olivier
 */
class GenericGF {
	/**
	 * Create a representation of GF(size) using the given primitive polynomial.
	 *
	 * @param primitive irreducible polynomial whose coefficients are represented by
	 *  the bits of an int, where the least-significant bit represents the constant
	 *  coefficient
	 * @param size the size of the field
	 * @param b the factor b in the generator polynomial can be 0- or 1-based
	 *  (g(x) = (x+a^b)(x+a^(b+1))...(x+a^(b+2t-1))).
	 *  In most cases it should be 1, but for QR code it is 0.
	 */
	constructor(primitive, size, b) {
		this.primitive = primitive;
		this.size = size;
		this.generatorBase = b;

		this.expTable = new Int32Array(size);
		this.logTable = new Int32Array(size);
		let x = 1;
		for (let i = 0; i < size; i++) {
			this.expTable[i] = x;
			x *= 2; // we're assuming the generator alpha is 2
			if (x >= size) {
				x ^= primitive;
				x &= size - 1;
			}
		}
		for (let i = 0; i < size - 1; i++) {
			this.logTable[this.expTable[i]] = i;
		}
		// logTable[0] == 0 but this should never be used
		this.zero = new GenericGFPoly(this, new Int32Array([0]));
		this.one = new GenericGFPoly(this, new Int32Array([1]));
	}

	getZero() {
		return this.zero;
	}

	getOne() {
		return this.one;
	}

	/**
	 * @return the monomial representing coefficient * x^degree
	 */
	buildMonomial(degree, coefficient) {
		if (degree < 0) {
			throw new ReedSolomonException();
		}
		if (coefficient === 0) {
			return this.zero;
		}
		const coefficients = new Int32Array(degree + 1);
		coefficients[0] = coefficient;
		return new GenericGFPoly(this, coefficients);
	}

	/**
	 * Implements both addition and subtraction -- they are the same in GF(size).
	 *
	 * @return sum/difference of a and b
	 */
	static addOrSubtract(a, b) {
		return a ^ b;
	}

	/**
	 * @return 2 to the power of a in GF(size)
	 */
	exp(a) {
		return this.expTable[a];
	}

	/**
	 * @return base 2 log of a in GF(size)
	 */
	log(a) {
		if (a === 0) {
			throw new ReedSolomonException();
		}
		// console.log(this.logTable);
		return this.logTable[a];
	}

	/**
	 * @return multiplicative inverse of a
	 */
	inverse(a) {
		if (a === 0) {
			throw new ReedSolomonException();
		}
		return this.expTable[this.size - this.logTable[a] - 1];
	}

	/**
	 * @return product of a and b in GF(size)
	 */
	multiply(a, b) {
		if (a === 0 || b === 0) {
			return 0;
		}
		return this.expTable[(this.logTable[a] + this.logTable[b]) % (this.size - 1)];
	}

	getSize() {
		return this.size;
	}

	getGeneratorBase() {
		return this.generatorBase;
	}

	toString() {
		return 'GF(0x' + this.primitive.toString(16) + ',' + this.size + ')';
	}
}

/**
 * <p>Represents a polynomial whose coefficients are elements of a GF.
 * Instances of this class are immutable.</p>
 *
 * <p>Much credit is due to William Rucklidge since portions of this code are an indirect
 * port of his C++ Reed-Solomon implementation.</p>
 *
 * @author Sean Owen
 */
class GenericGFPoly {
	/**
	 * @param field the {@link GenericGF} instance representing the field to use
	 * to perform computations
	 * @param coefficients coefficients as ints representing elements of GF(size), arranged
	 * from most significant (highest-power term) coefficient to least significant
	 * @throws IllegalArgumentException if argument is null or empty,
	 * or if leading coefficient is 0 and this is not a
	 * constant polynomial (that is, it is not the monomial "0")
	 */
	constructor(field, coefficients) {
		if (coefficients.length === 0 || !field) {
			throw new ReedSolomonException();
		}
		this.field = field;
		const coefficientsLength = coefficients.length;
		if (coefficientsLength > 1 && coefficients[0] === 0) {
			// Leading term must be non-zero for anything except the constant polynomial "0"
			let firstNonZero = 1;
			while (firstNonZero < coefficientsLength && coefficients[firstNonZero] === 0) {
				firstNonZero++;
			}
			if (firstNonZero === coefficientsLength) {
				this.coefficients = new Int32Array([0]);
			} else {
				const len = coefficientsLength - firstNonZero;
				const newCoefficients = new Int32Array(len);
				for (let i = 0; i < len; i++) {
					newCoefficients[i] = coefficients[i + firstNonZero];
				}
				newCoefficients.set(coefficients.slice(firstNonZero), 0);
				// System.arraycopy(coefficients,
				//     firstNonZero,
				//     this.coefficients,
				//     0,
				//     this.coefficients.length);
				this.coefficients = newCoefficients;
			}
		} else {
			this.coefficients = coefficients;
		}
	}

	getCoefficients() {
		return this.coefficients;
	}

	/**
	 * @return degree of this polynomial
	 */
	getDegree() {
		return this.coefficients.length - 1;
	}

	/**
	 * @return true iff this polynomial is the monomial "0"
	 */
	isZero() {
		return this.coefficients[0] === 0;
	}

	/**
	 * @return coefficient of x^degree term in this polynomial
	 */
	getCoefficient(degree) {
		// console.log('getCoefficient this.coefficients degree:' + degree);
		// console.log(this.coefficients);
		return this.coefficients[this.coefficients.length - 1 - degree];
	}

	/**
	 * @return evaluation of this polynomial at a given point
	 */
	evaluateAt(a) {
		const field = this.field;
		if (a === 0) {
			// Just return the x^0 coefficient
			return this.getCoefficient(0);
		}
		if (a === 1) {
			// Just the sum of the coefficients
			let result = 0;
			for (const coefficient of this.coefficients) {
				result = GenericGF.addOrSubtract(result, coefficient);
			}
			return result;
		}
		let result = this.coefficients[0];
		for (const coefficient of this.coefficients) {
			result = GenericGF.addOrSubtract(field.multiply(a, result), coefficient);
		}
		return result;
	}

	addOrSubtract(other) {
		const field = this.field;
		if (field !== other.field) {
			throw new ReedSolomonException('GenericGFPolys do not have same GenericGF field');
		}
		if (this.isZero()) {
			return other;
		}
		if (other.isZero()) {
			return this;
		}
		const isSmallLonger = this.coefficients.length > other.coefficients.length;
		// const coefficients = this.coefficients;
		const smallerCoefficients = isSmallLonger ? other.coefficients : this.coefficients;
		const largerCoefficients = isSmallLonger ? this.coefficients : other.coefficients;
		// if (smallerCoefficients.length > largerCoefficients.length) {
		// 	const temp = smallerCoefficients;
		// 	smallerCoefficients = largerCoefficients;
		// 	largerCoefficients = temp;
		// }
		const lenLarge = largerCoefficients.length;
		const sumDiff = new Int32Array(lenLarge);
		const lengthDiff = lenLarge - smallerCoefficients.length;
		// Copy high-order terms only found in higher-degree polynomial's coefficients
		sumDiff.set(largerCoefficients.slice(0, lengthDiff), 0);
		// System.arraycopy(largerCoefficients, 0, sumDiff, 0, lengthDiff);
		// console.log('lengthDiff:' + lengthDiff + '/sumDiff:' + JSON.stringify(sumDiff));
		for (let i = lengthDiff; i < lenLarge; i++) {
			sumDiff[i] = GenericGF.addOrSubtract(smallerCoefficients[i - lengthDiff], largerCoefficients[i]);
		}
		// console.log(sumDiff);
		return new GenericGFPoly(field, sumDiff);
	}

	multiply(param) {
		return typeof param === 'number' ? this.multiplyScalar(param) : this.multiplyGF(param);
	}
	multiplyGF(other) {
		const field = this.field;
		const aCoefficients = this.coefficients;
		if (field !== other.field) {
			throw new ReedSolomonException('GenericGFPolys do not have same GenericGF field');
		}
		if (this.isZero() || other.isZero()) {
			return field.getZero();
		}
		// const aCoefficients = coefficients;largerCoefficients
		const aLength = aCoefficients.length;
		const bCoefficients = other.coefficients;
		const bLength = bCoefficients.length;
		const product = new Int32Array(aLength + bLength - 1);
		for (let i = 0; i < aLength; i++) {
			const aCoeff = aCoefficients[i];
			for (let j = 0; j < bLength; j++) {
				product[i + j] = GenericGF.addOrSubtract(product[i + j], field.multiply(aCoeff, bCoefficients[j]));
			}
		}
		return new GenericGFPoly(field, product);
	}

	multiplyScalar(scalar) {
		const field = this.field;
		const coefficients = this.coefficients;
		if (scalar === 0) {
			return field.getZero();
		}
		if (scalar === 1) {
			return this;
		}
		const size = coefficients.length;
		const product = new Int32Array(size);
		for (let i = 0; i < size; i++) {
			product[i] = field.multiply(coefficients[i], scalar);
		}
		return new GenericGFPoly(field, product);
	}

	multiplyByMonomial(degree, coefficient) {
		const field = this.field;
		const coefficients = this.coefficients;
		if (degree < 0) {
			throw new ReedSolomonException();
		}
		if (coefficient === 0) {
			return this, field.getZero();
		}
		const size = coefficients.length;
		const product = new Int32Array(size + degree);
		for (let i = 0; i < size; i++) {
			product[i] = field.multiply(coefficients[i], coefficient);
		}
		return new GenericGFPoly(field, product);
	}

	divide(other) {
		const field = this.field;
		if (field !== other.field) {
			throw new ReedSolomonException('GenericGFPolys do not have same GenericGF field');
		}
		if (other.isZero()) {
			throw new ReedSolomonException('Divide by 0');
		}
		let quotient = field.getZero();
		let remainder = this;
		const denominatorLeadingTerm = other.getCoefficient(other.getDegree());
		const inverseDenominatorLeadingTerm = field.inverse(denominatorLeadingTerm);
		while (remainder.getDegree() >= other.getDegree() && !remainder.isZero()) {
			const degreeDifference = remainder.getDegree() - other.getDegree();
			const scale = field.multiply(remainder.getCoefficient(remainder.getDegree()), inverseDenominatorLeadingTerm);
			const term = other.multiplyByMonomial(degreeDifference, scale);
			const iterationQuotient = field.buildMonomial(degreeDifference, scale);
			quotient = quotient.addOrSubtract(iterationQuotient);
			remainder = remainder.addOrSubtract(term);
		}
		return [quotient, remainder];
	}

	toString() {
		const field = this.field;
		const coefficients = this.coefficients;
		if (this.isZero()) {
			return '0';
		}
		const result = []; //new StringBuilder(8 * getDegree());
		for (let degree = this.getDegree(); degree >= 0; degree--) {
			const coefficient = this.getCoefficient(degree);
			if (coefficient !== 0) {
				if (coefficient < 0) {
					if (degree === this.getDegree()) {
						result.push('-');
					} else {
						result.push(' - ');
					}
					coefficient = -coefficient;
				} else {
					if (result.length > 0) {
						result.push(' + ');
					}
				}
				if (degree === 0 || coefficient !== 1) {
					const alphaPower = field.log(coefficient);
					if (alphaPower === 0) {
						result.push('1');
					} else if (alphaPower === 1) {
						result.push('a');
					} else {
						result.push('a^');
						result.push(alphaPower);
					}
				}
				if (degree !== 0) {
					if (degree === 1) {
						result.push('x');
					} else {
						result.push('x^');
						result.push(degree);
					}
				}
			}
		}
		return result.join('');
	}
}

/**
 * <p>Implements Reed-Solomon decoding, as the name implies.</p>
 *
 * <p>The algorithm will not be explained here, but the following references were helpful
 * in creating this implementation:</p>
 *
 * <ul>
 * <li>Bruce Maggs.
 * <a href="http://www.cs.cmu.edu/afs/cs.cmu.edu/project/pscico-guyb/realworld/www/rs_decode.ps">
 * "Decoding Reed-Solomon Codes"</a> (see discussion of Forney's Formula)</li>
 * <li>J.I. Hall. <a href="www.mth.msu.edu/~jhall/classes/codenotes/GRS.pdf">
 * "Chapter 5. Generalized Reed-Solomon Codes"</a>
 * (see discussion of Euclidean algorithm)</li>
 * </ul>
 *
 * <p>Much credit is due to William Rucklidge since portions of this code are an indirect
 * port of his C++ Reed-Solomon implementation.</p>
 *
 * @author Sean Owen
 * @author William Rucklidge
 * @author sanfordsquires
 */
class ReedSolomonDecoder {
	constructor(field, isSloppy = true) {
		this.field = field;
		this.isSloppy = isSloppy;
	}

	/**
	 * <p>Decodes given set of received codewords, which include both data and error-correction
	 * codewords. Really, this means it uses Reed-Solomon to detect and correct errors, in-place,
	 * in the input.</p>
	 *
	 * @param received data and error-correction codewords
	 * @param twoS number of error-correction codewords available
	 * @throws ReedSolomonException if decoding fails for any reason
	 */
	decode(received, twoS) {
		const field = this.field;
		const poly = new GenericGFPoly(field, received);
		const syndromeCoefficients = new Int32Array(twoS);
		const scLen = syndromeCoefficients.length;
		let noError = true;
		const len = received.length;
		const generatorBase = field.getGeneratorBase();
		for (let i = 0; i < twoS; i++) {
			const evaled = poly.evaluateAt(field.exp(i + generatorBase));
			syndromeCoefficients[scLen - 1 - i] = evaled;
			if (evaled !== 0) {
				noError = false;
			}
		}
		// console.log('decode syndromeCoefficients noError:' + noError + '/twoS:' + twoS);
		// console.log(syndromeCoefficients);
		// console.log('decode received len:' + len);
		// console.log(received);
		if (noError) {
			return received.slice(0, len - twoS);
		}
		const syndrome = new GenericGFPoly(field, syndromeCoefficients);
		const [sigma, omega] = this.runEuclideanAlgorithm(field.buildMonomial(twoS, 1), syndrome, twoS);
		const errorLocations = this.findErrorLocations(sigma);
		// console.log('decode errorLocations sigma:' + sigma);
		// console.log(errorLocations);
		const errorMagnitudes = this.findErrorMagnitudes(omega, errorLocations);
		for (let i = 0; i < errorLocations.length; i++) {
			const errorLocation = errorLocations[i];
			const log = field.log(errorLocation);
			const position = len - 1 - log;
			// console.log('errorLocation:' + errorLocation + '/log:' + log + '/position:' + position);
			if (position < 0) {
				if (this.isSloppy) {
					continue;
				}
				throw new ReedSolomonException('Bad error location');
			}
			received[position] = GenericGF.addOrSubtract(received[position], errorMagnitudes[i]);
		}
		return received.slice(0, len - twoS);
	}

	runEuclideanAlgorithm(a, b, R) {
		const field = this.field;
		// Assume a's degree is >= b's
		if (a.getDegree() < b.getDegree()) {
			const temp = a;
			a = b;
			b = temp;
		}
		let rLast = a;
		let r = b;
		let tLast = field.getZero();
		let t = field.getOne();
		// Run Euclidean algorithm until r's degree is less than R/2
		while (r.getDegree() >= R / 2) {
			const rLastLast = rLast;
			const tLastLast = tLast;
			rLast = r;
			tLast = t;
			// Divide rLastLast by rLast, with quotient in q and remainder in r
			if (rLast.isZero()) {
				// Oops, Euclidean algorithm already terminated?
				throw new ReedSolomonException('r_{i-1} was zero');
			}
			r = rLastLast;
			let q = field.getZero();
			const rLastDegree = rLast.getDegree();
			const denominatorLeadingTerm = rLast.getCoefficient(rLastDegree);
			const dltInverse = field.inverse(denominatorLeadingTerm);
			while (r.getDegree() >= rLastDegree && !r.isZero()) {
				const rDegree = r.getDegree();
				const degreeDiff = rDegree - rLastDegree;
				const scale = field.multiply(r.getCoefficient(rDegree), dltInverse);
				q = q.addOrSubtract(field.buildMonomial(degreeDiff, scale));
				r = r.addOrSubtract(rLast.multiplyByMonomial(degreeDiff, scale));
			}
			t = q.multiply(tLast).addOrSubtract(tLastLast);
			if (r.getDegree() >= rLastDegree) {
				throw new ReedSolomonException('Division algorithm failed to reduce polynomial?');
			}
		}
		const sigmaTildeAtZero = t.getCoefficient(0);
		if (sigmaTildeAtZero === 0) {
			console.log('runEuclideanAlgorithm a b R t sigmaTildeAtZero:' + sigmaTildeAtZero);
			console.log(a);
			console.log(b);
			console.log(R);
			console.log(t);
			throw new ReedSolomonException('sigmaTilde(0) was zero');
		}
		const inverse = field.inverse(sigmaTildeAtZero);
		const sigma = t.multiply(inverse);
		const omega = r.multiply(inverse);
		return [sigma, omega];
	}

	findErrorLocations(errorLocator) {
		const field = this.field;
		// This is a direct application of Chien's search
		const numErrors = errorLocator.getDegree();
		if (numErrors === 1) {
			// shortcut
			// console.log('findErrorLocations errorLocator.getCoefficient(1):' + errorLocator.getCoefficient(1) + '/numErrors:' + numErrors);
			return new Int32Array([errorLocator.getCoefficient(1)]);
		}
		const result = new Int32Array(numErrors);
		let e = 0;
		const size = field.getSize();
		// console.log('findErrorLocations size:' + size + '/numErrors:' + numErrors);
		for (let i = 1; i < size && e < numErrors; i++) {
			if (errorLocator.evaluateAt(i) === 0) {
				result[e] = field.inverse(i);
				e++;
			}
		}
		if (e !== numErrors && !this.isSloppy) {
			throw new ReedSolomonException('Error locator degree does not match number of roots');
		}

		// console.log(result);
		return result.slice(0, e);
	}

	findErrorMagnitudes(errorEvaluator, errorLocations) {
		const field = this.field;
		// This is directly applying Forney's Formula
		const s = errorLocations.length;
		const result = new Int32Array(s);
		for (let i = 0; i < s; i++) {
			const xiInverse = field.inverse(errorLocations[i]);
			let denominator = 1;
			for (let j = 0; j < s; j++) {
				if (i !== j) {
					//denominator = field.multiply(denominator,
					//    GenericGF.addOrSubtract(1, field.multiply(errorLocations[j], xiInverse)));
					// Above should work but fails on some Apple and Linux JDKs due to a Hotspot bug.
					// Below is a funny-looking workaround from Steven Parkes
					const term = field.multiply(errorLocations[j], xiInverse);
					const termPlus1 = (term & 0x1) === 0 ? term | 1 : term & ~1;
					denominator = field.multiply(denominator, termPlus1);
				}
			}
			result[i] = field.multiply(errorEvaluator.evaluateAt(xiInverse), field.inverse(denominator));
			if (field.getGeneratorBase() !== 0) {
				result[i] = field.multiply(result[i], xiInverse);
			}
		}
		return result;
	}
}
/**
 * <p>Implements Reed-Solomon encoding, as the name implies.</p>
 *
 * @author Sean Owen
 * @author William Rucklidge
 */
class ReedSolomonEncoder {
	constructor(field) {
		this.field = field;
		this.cachedGenerators = [];
		const coefficients = new Int32Array([1]);
		this.cachedGenerators.push(new GenericGFPoly(field, coefficients));
	}

	buildGenerator(degree) {
		const field = this.field;
		const cachedGenerators = this.cachedGenerators;
		const len = cachedGenerators.length;
		if (degree >= len) {
			const lastIndex = len - 1;
			const generatorBase = field.getGeneratorBase();
			let lastGenerator = cachedGenerators[lastIndex];
			for (let d = lastIndex; d < degree; d++) {
				const coefficients = new Int32Array([1, field.exp(d + generatorBase)]);
				const nextGenerator = lastGenerator.multiply(new GenericGFPoly(field, coefficients));
				cachedGenerators.push(nextGenerator);
				lastGenerator = nextGenerator;
			}
		}
		return cachedGenerators[degree];
	}

	encode(toEncode, ecBytes) {
		const field = this.field;
		if (ecBytes === 0) {
			throw new ReedSolomonException('No error correction bytes');
		}
		const size = toEncode.length;
		let dataBytes = size - ecBytes;
		if (dataBytes <= 0) {
			// console.log('size:' + size + '/ecBytes:' + ecBytes);
			throw new ReedSolomonException('No data bytes provided');
		}
		const generator = this.buildGenerator(ecBytes);
		const infoCoefficients = new Int32Array(dataBytes);
		infoCoefficients.set(toEncode.slice(0, dataBytes), 0);
		// System.arraycopy(toEncode, 0, infoCoefficients, 0, dataBytes);
		const info = new GenericGFPoly(field, infoCoefficients);
		const info2 = info.multiplyByMonomial(ecBytes, 1);
		const remainder = info2.divide(generator)[1];
		const coefficients = remainder.getCoefficients();
		const numZeroCoefficients = ecBytes - coefficients.length;
		const result = new Int32Array(size);
		result.fill(0);
		result.set(infoCoefficients, 0);
		const offset = dataBytes + numZeroCoefficients;
		result.set(coefficients, offset);
		// console.log('encode toEncode');
		// console.log(toEncode);
		// console.log('encode result');
		// console.log(result);
		// System.arraycopy(coefficients, 0, toEncode, dataBytes + numZeroCoefficients, coefficients.length);
		return result;
	}
}
/**
 * <p>Thrown when an exception occurs during Reed-Solomon decoding, such as when
 * there are too many errors to correct.</p>
 *
 * @author Sean Owen
 */
export class ReedSolomonException extends Error {
	constructor(message) {
		super(message);
	}
}

export const AZTEC_DATA_12 = { primitive: 0x1069, bitNum: 12, b: 1 }; // x^12 + x^6 + x^5 + x^3 + 1 = 1 0000 0110 1001=4201=0x1069 r=12
export const AZTEC_DATA_10 = { primitive: 0x409, bitNum: 10, b: 1 }; // x^10 + x^3 + 1 =100 0000 1001=1033=0x409 r=10
export const AZTEC_DATA_6 = { primitive: 0x43, bitNum: 6, b: 1 }; // x^6 + x + 1=100 0011=67=0x43
export const AZTEC_PARAM = { primitive: 0x13, bitNum: 4, b: 1 }; // x^4 + x + 1=1 0011 = 19=0x13
export const QR_CODE_FIELD_256 = { primitive: 0x011d, bitNum: 8, b: 0 }; // x^8 + x^4 + x^3 + x^2 + 1=1 0001 1101=281=0x011d
export const DATA_MATRIX_FIELD_256 = { primitive: 0x012d, bitNum: 8, b: 1 }; // x^8 + x^5 + x^3 + x^2 + 1=1 0010 1101=297=0x012d
export const AZTEC_DATA_8 = DATA_MATRIX_FIELD_256;
export const MAXICODE_FIELD_64 = AZTEC_DATA_6;
const ByteAs16bit = { primitive: 0x10029, bitNum: 16, b: 1 }; // D^16+D^5+D^3+D^2+1= 1 0000 0000 0010 1101 =65581=0x10029 r=16;
const ByteAs14bit = { primitive: 0x402b, bitNum: 14, b: 1 }; // D^14+D^5+D^3+D+1= 100 0000 0010 1011 =16427=0x402b r=14;
const ByteAs12bit = AZTEC_DATA_12;
const ByteAs10bit = AZTEC_DATA_10;
const ByteAs8bit = AZTEC_DATA_8;
const ByteAs6bit = AZTEC_DATA_6;
const ByteAs4bit = AZTEC_PARAM;

export const presets = {
	ByteAs16bit,
	ByteAs14bit,
	ByteAs12bit,
	ByteAs10bit,
	ByteAs8bit,
	ByteAs6bit,
	ByteAs4bit,
	AZTEC_DATA_12,
	AZTEC_DATA_10,
	AZTEC_DATA_6,
	AZTEC_PARAM,
	QR_CODE_FIELD_256,
	DATA_MATRIX_FIELD_256,
	AZTEC_DATA_8,
	MAXICODE_FIELD_64,
};
const Cache = {};
export class ReedSolomonES {
	static encodeRaw(i32a, errorCorrectionRedundantUnitCount, primitive, bitNum, b) {
		const key = JSON.stringify([primitive, bitNum, b, 'e']);
		let rsEncoder = Cache[key];
		if (!rsEncoder) {
			const pow = Math.pow(2, bitNum);
			const ggf = new GenericGF(primitive, pow, b);
			rsEncoder = new ReedSolomonEncoder(ggf);
			Cache[key] = rsEncoder;
		}
		return rsEncoder.encode(i32a, errorCorrectionRedundantUnitCount);
	}
	static decodeRaw(i32a, errorCorrectionRedundantUnitCount, primitive, bitNum, b, isSloppy) {
		const key = JSON.stringify([primitive, bitNum, b, 'd']);
		let rsDecoder = Cache[key];
		if (!rsDecoder) {
			const pow = Math.pow(2, bitNum);
			const ggf = new GenericGF(primitive, pow, b);
			rsDecoder = new ReedSolomonDecoder(ggf, isSloppy);
			Cache[key] = rsDecoder;
		}
		return rsDecoder.decode(i32a, errorCorrectionRedundantUnitCount);
	}
	static copyToU8a(i32a, bitNum) {
		const dataLen = i32a.length;
		let lastOne = '';
		const newDataLength = Math.ceil((dataLen * bitNum) / 8);
		// const has = (dataLen * 8) % bitNum;
		const u8a = new Uint8Array(newDataLength);
		let charCounter = 0;
		const fill = new Array(bitNum);
		fill.fill('0');
		const fill0 = fill.join('');
		// const a = [];
		// const b = [];
		for (let octet of i32a) {
			const bits = (fill0 + octet.toString(2)).slice(-1 * bitNum);
			const current = lastOne + bits;
			// b.push(current);

			const loopCount = Math.floor(current.length / 8);
			for (let i = 0; i < loopCount; i++) {
				const bitStr = current.substring(i * 8, (i + 1) * 8);
				u8a[charCounter] = parseInt(bitStr, 2);
				// a.push(bitStr);
				charCounter++;
			}
			lastOne = current.substring(loopCount * 8);
		}
		// console.log(b);
		// console.log(a);
		if (!lastOne) {
			const a = new Array(8 - lastOne.length);
			a.fill('0');
			u8a[charCounter] = parseInt(lastOne + a.join(''), 2);
		}
		return u8a;
	}
	static copyToI32a(u8a, bitNum, isFillBit = true) {
		const dataLen = u8a.length;
		let lastOne = '';
		let lcm = bitNum;
		const bitnumHarf = bitNum / 2;
		let add = 0;
		if (isFillBit && bitNum !== 8) {
			const mod = dataLen % bitnumHarf;
			add = mod > 0 ? bitnumHarf - mod : mod;
		}
		// console.log('lcm:' + lcm + '/add:' + add + '/' + bitnumHarf);
		const newDataLength = Math.ceil(((dataLen + add) * 8) / bitNum);
		const i32a = new Int32Array(newDataLength);
		const fill = new Array(8);
		fill.fill('0');
		const fill0 = fill.join('');
		let charCounter = 0;
		for (let octet of u8a) {
			const bits = (fill0 + octet.toString(2)).slice(-8);
			const current = lastOne + bits;
			const loopCount = Math.floor(current.length / bitNum);
			for (let i = 0; i < loopCount; i++) {
				const bitStr = current.substring(i * bitNum, (i + 1) * bitNum);
				i32a[charCounter] = parseInt(bitStr, 2);
				charCounter++;
			}
			lastOne = current.substring(loopCount * bitNum);
		}
		// console.log('lastOne:' + lastOne + '/newDataLength:' + newDataLength + '/bitNum:' + bitNum + '/dataLen:' + dataLen);
		if (lastOne) {
			const a = new Array(bitNum - lastOne.length);
			a.fill('0');
			const bits = lastOne + a.join('');
			// console.log('bits:' + bits);
			i32a[charCounter] = parseInt(bits, 2);
		}
		// console.log(i32a);
		return i32a;
	}
	//2^{r}>N>K>0 ,t=(N-K)/2
	static encode(u8a, presetName, errorCorrectionRetio) {
		const preset = presets[presetName];
		if (!preset) {
			return;
		}
		const bitNum = preset.bitNum;
		const primitive = preset.primitive;
		const b = preset.b;
		const i32a = ReedSolomonES.copyToI32a(u8a, bitNum, true);
		const newDataLength = i32a.length;

		const maxWordLength = Math.pow(2, bitNum);
		const retio = errorCorrectionRetio * 2 + 1;
		const N = Math.ceil(retio * newDataLength);
		const wordCount = Math.ceil(N / maxWordLength);

		const wordKtmp = Math.floor(maxWordLength / retio);
		const wordK = Math.ceil((wordKtmp > newDataLength ? wordKtmp : newDataLength) / wordCount);
		const result = new Int32Array(newDataLength * retio);
		let leftLength = newDataLength;
		// console.log(newDataLength + '/' + errorCorrectionRetio + '/retio:' + retio + '/wordK:' + wordK + '/wordCount:' + wordCount);
		for (let i = 0; i < wordCount; i++) {
			const wordKCurrent = leftLength > wordK ? wordK : leftLength;
			const errorCorrectionRedundantUnitCount = Math.floor(wordKCurrent * errorCorrectionRetio * 2);
			const tempN = wordKCurrent + errorCorrectionRedundantUnitCount;
			leftLength -= wordK;
			const start = i * wordK;
			if (start > newDataLength || wordKCurrent < retio) {
				break;
			}
			const end = start + wordKCurrent > newDataLength ? newDataLength : start + wordKCurrent;
			const toEncode = new Int32Array(tempN);
			toEncode.fill(0);
			toEncode.set(i32a.slice(start, end), 0);
			// console.log('toEncode');
			// console.log(toEncode);
			const encorded = ReedSolomonES.encodeRaw(toEncode, errorCorrectionRedundantUnitCount, primitive, bitNum, b);
			const offset = tempN * i;
			result.set(encorded, offset);
			// console.log('encorded offset:' + offset);
			// console.log(encorded);
			// console.log(result);
		}
		// console.log('encode presetName:' + presetName + '/preset:' + JSON.stringify(preset));
		return ReedSolomonES.copyToU8a(result, bitNum);
	}
	static decode(u8a, presetName, errorCorrectionRetio, isSloppy) {
		const preset = presets[presetName];
		if (!preset) {
			return;
		}
		const bitNum = preset.bitNum;
		const primitive = preset.primitive;
		const b = preset.b;
		const i32a = ReedSolomonES.copyToI32a(u8a, bitNum);
		const newDataLength = i32a.length;

		const maxWordLength = Math.pow(2, bitNum);
		const retio = errorCorrectionRetio * 2 + 1;
		const K = Math.ceil(newDataLength / retio);
		const wordCount = Math.ceil(newDataLength / maxWordLength);

		const wordKtmp = Math.floor(maxWordLength / retio);
		const wordK = Math.ceil((wordKtmp > K ? wordKtmp : K) / wordCount);
		const result = new Int32Array(K);
		let leftLength = K;
		for (let i = 0; i < wordCount; i++) {
			const wordKCurrent = leftLength > wordK ? wordK : leftLength;
			const errorCorrectionRedundantUnitCount = Math.floor(wordKCurrent * errorCorrectionRetio * 2);
			const tempN = wordKCurrent + errorCorrectionRedundantUnitCount;
			leftLength -= wordKCurrent;
			const start = i * tempN;
			const end = start + tempN > newDataLength ? newDataLength : start + tempN;
			const toDecode = new Int32Array(tempN);
			const na = i32a.slice(start, end);
			console.log('i:' + i + '/na:' + na.length + '/tempN:' + tempN);
			if (na.length < 1) {
				break;
			}
			toDecode.set(na, 0);
			const decorded = ReedSolomonES.decodeRaw(toDecode, errorCorrectionRedundantUnitCount, primitive, bitNum, b, isSloppy);
			const offset = wordK * i;
			result.set(decorded, offset);
		}
		console.log('decode presetName:' + presetName + '/preset:' + JSON.stringify(preset));
		return ReedSolomonES.copyToU8a(result, bitNum);
	}
}
