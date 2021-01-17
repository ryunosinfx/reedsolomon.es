const max = 50000;
const RB64Regex = /^[0-9a-zA-Z/\+]+[=]{0,3}$/;
const td = new TextDecoder('utf-8');
const te = new TextEncoder();
export class B64Util {
	static from64(d) {
		const a = new Uint8Array(1);
		a.fill(0);
		const b = B64Util.b64ToU8a(d);
		const u8a = b.length % 2 ? B64Util.joinU8as([b, a]) : b;
		const u16a = new Uint16Array(u8a.buffer);
		const l = u16a.length;
		const c = Math.ceil(l / max);
		const r = [];
		for (let j = 0; j < c; j++) {
			const start = max * j;
			const size = l - start;
			const p = size > max ? max : size > 0 ? size : l;
			const u = u16a.slice(start, start + p);
			r.push(String.fromCharCode(...u));
		}
		return r.join('');
	}
	static to64(s) {
		const len = s.length;
		const pageNum = Math.ceil(len / max);
		const results = [];
		for (let j = 0; j < pageNum; j++) {
			const start = max * j;
			const size = len - start;
			const p = size > max ? max : size > 0 ? size : len;
			const end = start + p;
			const input = s.substring(start, end);
			// const u = new Uint16Array(p);
			// for (let i = 0; i < p; i++) {
			// 	u[i] = input.charCodeAt(i);
			// }
			const ab = te.encode(input);
			const c = String.fromCharCode(...new Uint8Array(ab));
			results.push(c);
		}
		return btoa(results.join(''));
	}
	static b64ToU8a(d) {
		const a = atob(d);
		const b = new Uint8Array(a.length);
		for (let i = 0; i < b.length; i++) {
			try {
				b[i] = a.charCodeAt(i);
			} catch (e) {
				console.log(i);
				console.log(e);
			}
		}
		return b;
	}
	static u8a2b64(u8a) {
		const bs = u8a ? B64Util.u8a2bs(u8a) : null;
		return bs ? btoa(bs) : null;
	}
	static s2u8a(s) {
		const d = B64Util.to64(s);
		return B64Util.b64ToU8a(d);
	}
	static s2hex(s) {
		const d = B64Util.to64(s);
		const hex = B64Util.b64toHex(d);
		return hex;
	}
	static hex2s(hex) {
		const u8a = B64Util.hex2u8a(hex);
		const d = B64Util.aToB64(u8a.buffer);
		return B64Util.from64(d);
	}
	static b64uToAb(b) {
		const d = B64Util.toB64(b);
		return B64Util.b64ToU8a(d).buffer;
	}
	static u8a2bs(u8a) {
		const r = [];
		for (let e of u8a) {
			r.push(String.fromCharCode(e));
		}
		return r.join('');
	}
	static hex2u8a(hex) {
		return new Uint8Array(
			hex.match(/[0-9a-f]{2}/gi).map((h) => {
				return parseInt(h, 16);
			})
		);
	}

	static ab2bs(ab) {
		return B64Util.u8a2bs(new Uint8Array(ab));
	}
	static aToB64(ai) {
		const ab = ai.buffer ? ai.buffer : ai;
		return btoa(B64Util.ab2bs(ab));
	}
	static aToB64u(ai) {
		const b = B64Util.aToB64(ai);
		return B64Util.toB64u(b);
	}
	static b64toHex(b64) {
		const u8a = B64Util.b64ToU8a(b64);
		return B64Util.aToHex(u8a);
	}
	static aToHex(ai) {
		const u8a = ai.buffer ? new Uint8Array(ai.buffer) : new Uint8Array(ai);
		const rl = [];
		for (let i of u8a) {
			const a = i.toString(16);
			rl.push(('00' + a).slice(-2));
		}
		return rl.join('');
	}
	static bs2u8a(bs) {
		const l = bs.length;
		const a = new Uint8Array(new ArrayBuffer(l));
		for (let i = 0; i < l; i++) {
			a[i] = bs.charCodeAt(i);
		}
		return a;
	}
	static isB64(d) {
		return d && typeof d === 'string' && d.length % 4 === 0 && RB64Regex.test(d);
	}
	static u8aToUtf8(u8a) {
		return td.decode(u8a.buffer);
	}
	static bs2utf8(bs) {
		const u8a = B64Util.bs2u8a(bs);
		return td.decode(u8a.buffer);
	}
	static dataURI2bs(dURI) {
		return atob(dURI.split(',')[1]);
	}
	static dataURI2u8a(dURI) {
		return B64Util.bs2u8a(atob(dURI.split(',')[1]));
	}
	static ab2dataURI(ab, type = 'application/octet-stream') {
		const b = btoa(B64Util.ab2bs(ab));
		return 'data:' + type + ';base64,' + b;
	}
	static joinU8as(u8as) {
		let l = 0;
		for (let u8a of u8as) {
			l += u8a.length;
		}
		const r = new Uint8Array(l);
		let s = 0;
		for (let u8a of u8as) {
			r.set(u8a, s);
			s += u8a.length;
		}
		return r;
	}
	static toB64u(b) {
		return b ? b.split('+').join('-').split('/').join('_').split('=').join('') : b;
	}
	static toB64(b64u) {
		const l = b64u.length;
		const c = l % 4 > 0 ? 4 - (l % 4) : 0;
		let b = b64u.split('-').join('+').split('_').join('/');

		for (let i = 0; i < c; i++) {
			b += '=';
		}
		return b;
	}
	static async sig(u8a) {
		const bs = B64Util.u8a2bs(u8a);
		return Hasher.sha256(bs, 1, 'hex');
	}
}
