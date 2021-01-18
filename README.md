# reedsolomon.es
This is Reed-Solomon coding module which is porting from [Zxing](https://github.com/zxing/zxing) open source project written in Java. 


# Live demo

https://ryunosinfx.github.io/reedsolomon.es/index.html

# usage

## encode 
```EJS
<script type="module" src="./ReedSolomon.js">
  // t=(N-K)/2 (ByteAs8bit:max 0.97)
  const errorCrrectionReduntantRetio = 0.97;
  
  // ByteAs4bit*,ByteAs6bit*,ByteAs8bit,ByteAs10bit,ByteAs12bit,ByteAs14bit,ByteAs16bit*,QR_CODE_FIELD_256
  // *:has error,unusualbe.
  const presetName = 'ByteAs8bit';
  
  // K = planeUint8array.length
  const planeUint8array = new Uint8Array([xxxxxxxxxxxxxxxxxxxxxxxxxx]);
  
  const encordedUint8Array 
    = ReedSolomonES.encode(planeUint8array, presetName, errorCrrectionReduntantRetio);
      
</script>
```
## decode 
```EJS
<script type="module" src="./ReedSolomon.js">
  // t=(N-K)/2 (ByteAs8bit:max 0.97).
  // It is must same setting on encording.
  const errorCrrectionReduntantRetio = 0.97;
  
  // ByteAs4bit*,ByteAs6bit*,ByteAs8bit,ByteAs10bit,ByteAs12bit,ByteAs14bit,ByteAs16bit*,QR_CODE_FIELD_256
  // *:has error,unusualbe. It is must same setting on encording.
  const presetName = 'ByteAs8bit';
  
  // N = encordedUint8Array.length
  const encordedUint8Array = new Uint8Array([yyyyyyyyyyyyyyyyyyyyyyy]); 
   
  // If you use ByteAs10bit,ByteAs12bit,ByteAs14bit. Then add [0] byte add satisfaction bit num by LCM;
  
  // Best effort mode:if errors are correctable,they are corrected.
  // But errors are not correctable,they are outputed.
  const decodedUint8Array 
    = ReedSolomonES.decode(encordedUint8Array, presetName, errorCrrectionReduntantRetio);
  
  // Strict mode:if errors are correctable,they are corrected.
  // But errors are not correctable,throw ReedSolomonException.
  const correctedUint8Array 
    = ReedSolomonES.decodeStrict(encordedUint8Array, presetName, errorCrrectionReduntantRetio);

</script>
```

## refarence（java）
Source:

https://github.com/zxing/zxing/tree/master/core/src/main/java/com/google/zxing/common/reedsolomon

JavaDoc:

https://zxing.github.io/zxing/apidocs/com/google/zxing/common/reedsolomon/package-frame.html


