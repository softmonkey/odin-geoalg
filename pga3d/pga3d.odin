// 3D Projective Geometric Algebra
// Based on code from https://bivector.net/tools.html?p=3&q=0&r=1

package bivector_pga3d

import "core:math"
import "core:strings"
import "core:fmt"
import "core:reflect"
import "core:io"
import "core:testing"

// 3DPGA type
pga3d :: distinct [16]f32

// Proc groups
add       :: proc { add_pga3d, add_multi_pga3d, adds_pga3d, sadd_pga3d, }
sub       :: proc { sub_pga3d, sub_multi_pga3d, ssub_pga3d, subs_pga3d, }
outer     :: proc { outer_pga3d, outer_multi_pga3d, }
mul       :: proc { mul_pga3d, mul_multi_pga3d, smul_pga3d, muls_pga3d, }
reverse   :: proc { reverse_pga3d, }
dual      :: proc { dual_pga3d, }
conj      :: proc { conjugate_pga3d, }
involute  :: proc { involute_pga3d, }
regress   :: proc { regress_pga3d, }
inner     :: proc { inner_pga3d, }
dot       :: proc { inner_pga3d, }

// fmt support
COLORED_FMT := true

@(private="file")
pga3d_basis_labels : [16]string = { "","e0","e1","e2","e3","e0e1","e0e2","e0e3","e1e2","e3e1","e2e3","e0e2e1","e0e1e3","e0e3e2","e1e2e3","e0e1e2e3" }

@(init, private="file")
init :: proc() {
  fmt.set_user_formatters(&pga3d_formatters)
  err := fmt.register_user_formatter(type_info_of(pga3d).id, User_Formatter)
  assert(err == .None)
}

@(private="file")
pga3d_formatters: map[typeid]fmt.User_Formatter

@(private="file")
User_Formatter :: proc(fi: ^fmt.Info, arg: any, verb: rune) -> bool {
  m := cast(^pga3d)arg.data
  switch verb {
    case 'v':
      io.write_string(fi.writer, fmt_pga3d(m))
    case 'f':
      fmt.fmt_array(fi, m, 16, size_of(f32), type_info_of(f32),'f')
    case:
      return false
    }
    return true
}

@(private="file")
fmt_pga3d :: proc(a: ^pga3d) -> string {
  builder, err := strings.builder_make_len_cap(0,256,context.temp_allocator)
  if err != nil { return "" }
  defer strings.builder_destroy(&builder)

  tmpl := COLORED_FMT ? "%s\x1B[38;5;81m%.3f\x1B[38;5;214m%s\033[0m" : "%s%.3f%s"

  for val,idx in a {
    if math.abs(val) >= math.F32_EPSILON {
      fmt.sbprintf( &builder, tmpl, 
        len(builder.buf)>0 ? ", " : "", 
        val, 
        pga3d_basis_labels[idx] )
    }
  }
  if len(builder.buf) == 0 { return "\x1B[38;5;81m0.000\033[0m" }
  return strings.to_string(builder)
}

// procs

// Reverse ~a
reverse_pga3d :: proc(a : pga3d) -> (res: pga3d) {
  res[0]   =  a[0]
  res[1]   =  a[1]
  res[2]   =  a[2]
  res[3]   =  a[3]
  res[4]   =  a[4]
  res[5]   = -a[5]
  res[6]   = -a[6]
  res[7]   = -a[7]
  res[8]   = -a[8]
  res[9]   = -a[9]
  res[10]  = -a[10]
  res[11]  = -a[11]
  res[12]  = -a[12]
  res[13]  = -a[13]
  res[14]  = -a[14]
  res[15]  =  a[15]
  return
}

// Dual !a
dual_pga3d :: proc(a : pga3d) -> (res: pga3d) {
  res[0]   =  a[15]
  res[1]   =  a[14]
  res[2]   =  a[13]
  res[3]   =  a[12]
  res[4]   =  a[11]
  res[5]   = -a[10]
  res[6]   = -a[9]
  res[7]   = -a[8]
  res[8]   = -a[7]
  res[9]   = -a[6]
  res[10]  = -a[5]
  res[11]  = -a[4]
  res[12]  = -a[3]
  res[13]  = -a[2]
  res[14]  = -a[1]
  res[15]  = -a[0]
  return
}

// Conjugate
conjugate_pga3d :: proc(a : pga3d) -> (res: pga3d) {
  res[0]   =  a[0]
  res[1]   = -a[1]
  res[2]   = -a[2]
  res[3]   = -a[3]
  res[4]   = -a[4]
  res[5]   = -a[5]
  res[6]   = -a[6]
  res[7]   = -a[7]
  res[8]   = -a[8]
  res[9]   = -a[9]
  res[10]  = -a[10]
  res[11]  =  a[11]
  res[12]  =  a[12]
  res[13]  =  a[13]
  res[14]  =  a[14]
  res[15]  =  a[15]
  return
}

// Involute
involute_pga3d :: proc(a: pga3d) -> (res: pga3d) {
  res[0]   =  a[0]
  res[1]   = -a[1]
  res[2]   = -a[2]
  res[3]   = -a[3]
  res[4]   = -a[4]
  res[5]   =  a[5]
  res[6]   =  a[6]
  res[7]   =  a[7]
  res[8]   =  a[8]
  res[9]   =  a[9]
  res[10]  =  a[10]
  res[11]  = -a[11]
  res[12]  = -a[12]
  res[13]  = -a[13]
  res[14]  = -a[14]
  res[15]  =  a[15]
  return
}

// Mul (a*b)
mul_pga3d :: proc(a: pga3d, b: pga3d) -> (res: pga3d) {
  res[0]=b[0]*a[0]+b[2]*a[2]+b[3]*a[3]+b[4]*a[4]-b[8]*a[8]-b[9]*a[9]-b[10]*a[10]-b[14]*a[14]
  res[1]=b[1]*a[0]+b[0]*a[1]-b[5]*a[2]-b[6]*a[3]-b[7]*a[4]+b[2]*a[5]+b[3]*a[6]+b[4]*a[7]+b[11]*a[8]+b[12]*a[9]+b[13]*a[10]+b[8]*a[11]+b[9]*a[12]+b[10]*a[13]+b[15]*a[14]-b[14]*a[15]
  res[2]=b[2]*a[0]+b[0]*a[2]-b[8]*a[3]+b[9]*a[4]+b[3]*a[8]-b[4]*a[9]-b[14]*a[10]-b[10]*a[14]
  res[3]=b[3]*a[0]+b[8]*a[2]+b[0]*a[3]-b[10]*a[4]-b[2]*a[8]-b[14]*a[9]+b[4]*a[10]-b[9]*a[14]
  res[4]=b[4]*a[0]-b[9]*a[2]+b[10]*a[3]+b[0]*a[4]-b[14]*a[8]+b[2]*a[9]-b[3]*a[10]-b[8]*a[14]
  res[5]=b[5]*a[0]+b[2]*a[1]-b[1]*a[2]-b[11]*a[3]+b[12]*a[4]+b[0]*a[5]-b[8]*a[6]+b[9]*a[7]+b[6]*a[8]-b[7]*a[9]-b[15]*a[10]-b[3]*a[11]+b[4]*a[12]+b[14]*a[13]-b[13]*a[14]-b[10]*a[15]
  res[6]=b[6]*a[0]+b[3]*a[1]+b[11]*a[2]-b[1]*a[3]-b[13]*a[4]+b[8]*a[5]+b[0]*a[6]-b[10]*a[7]-b[5]*a[8]-b[15]*a[9]+b[7]*a[10]+b[2]*a[11]+b[14]*a[12]-b[4]*a[13]-b[12]*a[14]-b[9]*a[15]
  res[7]=b[7]*a[0]+b[4]*a[1]-b[12]*a[2]+b[13]*a[3]-b[1]*a[4]-b[9]*a[5]+b[10]*a[6]+b[0]*a[7]-b[15]*a[8]+b[5]*a[9]-b[6]*a[10]+b[14]*a[11]-b[2]*a[12]+b[3]*a[13]-b[11]*a[14]-b[8]*a[15]
  res[8]=b[8]*a[0]+b[3]*a[2]-b[2]*a[3]+b[14]*a[4]+b[0]*a[8]+b[10]*a[9]-b[9]*a[10]+b[4]*a[14]
  res[9]=b[9]*a[0]-b[4]*a[2]+b[14]*a[3]+b[2]*a[4]-b[10]*a[8]+b[0]*a[9]+b[8]*a[10]+b[3]*a[14]
  res[10]=b[10]*a[0]+b[14]*a[2]+b[4]*a[3]-b[3]*a[4]+b[9]*a[8]-b[8]*a[9]+b[0]*a[10]+b[2]*a[14]
  res[11]=b[11]*a[0]-b[8]*a[1]+b[6]*a[2]-b[5]*a[3]+b[15]*a[4]-b[3]*a[5]+b[2]*a[6]-b[14]*a[7]-b[1]*a[8]+b[13]*a[9]-b[12]*a[10]+b[0]*a[11]+b[10]*a[12]-b[9]*a[13]+b[7]*a[14]-b[4]*a[15]
  res[12]=b[12]*a[0]-b[9]*a[1]-b[7]*a[2]+b[15]*a[3]+b[5]*a[4]+b[4]*a[5]-b[14]*a[6]-b[2]*a[7]-b[13]*a[8]-b[1]*a[9]+b[11]*a[10]-b[10]*a[11]+b[0]*a[12]+b[8]*a[13]+b[6]*a[14]-b[3]*a[15]
  res[13]=b[13]*a[0]-b[10]*a[1]+b[15]*a[2]+b[7]*a[3]-b[6]*a[4]-b[14]*a[5]-b[4]*a[6]+b[3]*a[7]+b[12]*a[8]-b[11]*a[9]-b[1]*a[10]+b[9]*a[11]-b[8]*a[12]+b[0]*a[13]+b[5]*a[14]-b[2]*a[15]
  res[14]=b[14]*a[0]+b[10]*a[2]+b[9]*a[3]+b[8]*a[4]+b[4]*a[8]+b[3]*a[9]+b[2]*a[10]+b[0]*a[14]
  res[15]=b[15]*a[0]+b[14]*a[1]+b[13]*a[2]+b[12]*a[3]+b[11]*a[4]+b[10]*a[5]+b[9]*a[6]+b[8]*a[7]+b[7]*a[8]+b[6]*a[9]+b[5]*a[10]-b[4]*a[11]-b[3]*a[12]-b[2]*a[13]-b[1]*a[14]+b[0]*a[15]
  return
}

mul_multi_pga3d :: proc(a: ..pga3d) -> (res: pga3d) {
  res = a[0]
  for i in a[1:] {
    res = mul_pga3d(res, i)
  }
  return
}

smul_pga3d :: proc(a: f32, b: pga3d) -> (res: pga3d) { return a*b }
muls_pga3d :: proc(a: pga3d, b: f32) -> (res: pga3d) { return a*b }

// Outer product (a ^ b)
outer_pga3d :: proc(a: pga3d, b: pga3d) -> (res: pga3d) {
  res[0]=b[0]*a[0]
  res[1]=b[1]*a[0]+b[0]*a[1]
  res[2]=b[2]*a[0]+b[0]*a[2]
  res[3]=b[3]*a[0]+b[0]*a[3]
  res[4]=b[4]*a[0]+b[0]*a[4]
  res[5]=b[5]*a[0]+b[2]*a[1]-b[1]*a[2]+b[0]*a[5]
  res[6]=b[6]*a[0]+b[3]*a[1]-b[1]*a[3]+b[0]*a[6]
  res[7]=b[7]*a[0]+b[4]*a[1]-b[1]*a[4]+b[0]*a[7]
  res[8]=b[8]*a[0]+b[3]*a[2]-b[2]*a[3]+b[0]*a[8]
  res[9]=b[9]*a[0]-b[4]*a[2]+b[2]*a[4]+b[0]*a[9]
  res[10]=b[10]*a[0]+b[4]*a[3]-b[3]*a[4]+b[0]*a[10]
  res[11]=b[11]*a[0]-b[8]*a[1]+b[6]*a[2]-b[5]*a[3]-b[3]*a[5]+b[2]*a[6]-b[1]*a[8]+b[0]*a[11]
  res[12]=b[12]*a[0]-b[9]*a[1]-b[7]*a[2]+b[5]*a[4]+b[4]*a[5]-b[2]*a[7]-b[1]*a[9]+b[0]*a[12]
  res[13]=b[13]*a[0]-b[10]*a[1]+b[7]*a[3]-b[6]*a[4]-b[4]*a[6]+b[3]*a[7]-b[1]*a[10]+b[0]*a[13]
  res[14]=b[14]*a[0]+b[10]*a[2]+b[9]*a[3]+b[8]*a[4]+b[4]*a[8]+b[3]*a[9]+b[2]*a[10]+b[0]*a[14]
  res[15]=b[15]*a[0]+b[14]*a[1]+b[13]*a[2]+b[12]*a[3]+b[11]*a[4]+b[10]*a[5]+b[9]*a[6]+b[8]*a[7]+b[7]*a[8]+b[6]*a[9]+b[5]*a[10]-b[4]*a[11]-b[3]*a[12]-b[2]*a[13]-b[1]*a[14]+b[0]*a[15]  
  return
}

outer_multi_pga3d :: proc(a: ..pga3d) -> (res: pga3d) {
  res = a[0]
  for i in a[1:] {
    res = outer_pga3d(res, i)
  }
  return
}

// Regressive product (a&b)
regress_pga3d :: proc(a: pga3d, b: pga3d) -> (res: pga3d) {
  res[15] = 1*(a[15]*b[15])
  res[14] = -1*(a[14]*-1*b[15]+a[15]*b[14]*-1)
  res[13] = -1*(a[13]*-1*b[15]+a[15]*b[13]*-1)
  res[12] = -1*(a[12]*-1*b[15]+a[15]*b[12]*-1)
  res[11] = -1*(a[11]*-1*b[15]+a[15]*b[11]*-1)
  res[10] = 1*(a[10]*b[15]+a[13]*-1*b[14]*-1-a[14]*-1*b[13]*-1+a[15]*b[10])
  res[9]  = 1*(a[9]*b[15]+a[12]*-1*b[14]*-1-a[14]*-1*b[12]*-1+a[15]*b[9])
  res[8]  = 1*(a[8]*b[15]+a[11]*-1*b[14]*-1-a[14]*-1*b[11]*-1+a[15]*b[8])
  res[7]  = 1*(a[7]*b[15]+a[12]*-1*b[13]*-1-a[13]*-1*b[12]*-1+a[15]*b[7])
  res[6]  = 1*(a[6]*b[15]-a[11]*-1*b[13]*-1+a[13]*-1*b[11]*-1+a[15]*b[6])
  res[5]  = 1*(a[5]*b[15]+a[11]*-1*b[12]*-1-a[12]*-1*b[11]*-1+a[15]*b[5])
  res[4]  = 1*(a[4]*b[15]-a[7]*b[14]*-1+a[9]*b[13]*-1-a[10]*b[12]*-1-a[12]*-1*b[10]+a[13]*-1*b[9]-a[14]*-1*b[7]+a[15]*b[4])
  res[3]  = 1*(a[3]*b[15]-a[6]*b[14]*-1-a[8]*b[13]*-1+a[10]*b[11]*-1+a[11]*-1*b[10]-a[13]*-1*b[8]-a[14]*-1*b[6]+a[15]*b[3])
  res[2]  = 1*(a[2]*b[15]-a[5]*b[14]*-1+a[8]*b[12]*-1-a[9]*b[11]*-1-a[11]*-1*b[9]+a[12]*-1*b[8]-a[14]*-1*b[5]+a[15]*b[2])
  res[1]  = 1*(a[1]*b[15]+a[5]*b[13]*-1+a[6]*b[12]*-1+a[7]*b[11]*-1+a[11]*-1*b[7]+a[12]*-1*b[6]+a[13]*-1*b[5]+a[15]*b[1])
  res[0]  = 1*(a[0]*b[15]+a[1]*b[14]*-1+a[2]*b[13]*-1+a[3]*b[12]*-1+a[4]*b[11]*-1+a[5]*b[10]+a[6]*b[9]+a[7]*b[8]+a[8]*b[7]+a[9]*b[6]+a[10]*b[5]-a[11]*-1*b[4]-a[12]*-1*b[3]-a[13]*-1*b[2]-a[14]*-1*b[1]+a[15]*b[0])
  return
}

// Inner product (a|b)
inner_pga3d :: proc(a: pga3d, b: pga3d) -> (res: pga3d) {
  res[0]=b[0]*a[0]+b[2]*a[2]+b[3]*a[3]+b[4]*a[4]-b[8]*a[8]-b[9]*a[9]-b[10]*a[10]-b[14]*a[14]
  res[1]=b[1]*a[0]+b[0]*a[1]-b[5]*a[2]-b[6]*a[3]-b[7]*a[4]+b[2]*a[5]+b[3]*a[6]+b[4]*a[7]+b[11]*a[8]+b[12]*a[9]+b[13]*a[10]+b[8]*a[11]+b[9]*a[12]+b[10]*a[13]+b[15]*a[14]-b[14]*a[15]
  res[2]=b[2]*a[0]+b[0]*a[2]-b[8]*a[3]+b[9]*a[4]+b[3]*a[8]-b[4]*a[9]-b[14]*a[10]-b[10]*a[14]
  res[3]=b[3]*a[0]+b[8]*a[2]+b[0]*a[3]-b[10]*a[4]-b[2]*a[8]-b[14]*a[9]+b[4]*a[10]-b[9]*a[14]
  res[4]=b[4]*a[0]-b[9]*a[2]+b[10]*a[3]+b[0]*a[4]-b[14]*a[8]+b[2]*a[9]-b[3]*a[10]-b[8]*a[14]
  res[5]=b[5]*a[0]-b[11]*a[3]+b[12]*a[4]+b[0]*a[5]-b[15]*a[10]-b[3]*a[11]+b[4]*a[12]-b[10]*a[15]
  res[6]=b[6]*a[0]+b[11]*a[2]-b[13]*a[4]+b[0]*a[6]-b[15]*a[9]+b[2]*a[11]-b[4]*a[13]-b[9]*a[15]
  res[7]=b[7]*a[0]-b[12]*a[2]+b[13]*a[3]+b[0]*a[7]-b[15]*a[8]-b[2]*a[12]+b[3]*a[13]-b[8]*a[15]
  res[8]=b[8]*a[0]+b[14]*a[4]+b[0]*a[8]+b[4]*a[14]
  res[9]=b[9]*a[0]+b[14]*a[3]+b[0]*a[9]+b[3]*a[14]
  res[10]=b[10]*a[0]+b[14]*a[2]+b[0]*a[10]+b[2]*a[14]
  res[11]=b[11]*a[0]+b[15]*a[4]+b[0]*a[11]-b[4]*a[15]
  res[12]=b[12]*a[0]+b[15]*a[3]+b[0]*a[12]-b[3]*a[15]
  res[13]=b[13]*a[0]+b[15]*a[2]+b[0]*a[13]-b[2]*a[15]
  res[14]=b[14]*a[0]+b[0]*a[14]
  res[15]=b[15]*a[0]+b[0]*a[15]
  return
}

inner_multi_pga3d :: proc(a: ..pga3d) -> (res: pga3d) {
  res = a[0]
  for i in a[1:] {
    res = inner_pga3d(res, i)
  }
  return
}

// Addition (a+b)
add_pga3d :: proc(a: pga3d, b: pga3d) -> (res: pga3d) { return a+b }

@(require_results) 
sadd_pga3d :: proc(a: f32, b: pga3d) -> (res: pga3d) {
  res = b 
  res[0] = a+b[0]
  return 
}

@(require_results) 
adds_pga3d :: proc(a: pga3d, b: f32) -> (res: pga3d) { 
  res = a
  res[0] = a[0]+b
  return
}

add_multi_pga3d :: proc(a: ..pga3d) -> (res: pga3d) {
  res = a[0]
  for i in a[1:] {
    res = add_pga3d(res, i)
  }
  return
}

// Subtraction (a-b)
@(require_results) sub_pga3d :: proc(a: pga3d, b: pga3d) -> (res: pga3d) { return a-b }

@(require_results) ssub_pga3d :: proc(a: f32, b: pga3d) -> (res: pga3d) { 
  res = -b
  res[0] = a-b[0]
  return
}

@(require_results) subs_pga3d :: proc(a: pga3d, b: f32) -> (res: pga3d) { 
  res = a
  res[0] = a[0]-b
  return
}

sub_multi_pga3d :: proc(a: ..pga3d) -> (res: pga3d) {
  res = a[0]
  for i in a[1:] {
    res = sub_pga3d(res, i)
  }
  return
}


norm :: proc(a: pga3d) -> f32 {
  return math.sqrt(math.abs(mul(a,conj(a))[0]))
}

inorm :: proc(a: pga3d) -> f32 {
  return norm( dual(a) )
}

normalized :: proc(a: pga3d) -> pga3d {
  return mul(a, 1.0/norm(a))
}

lerp :: proc(a: pga3d, b: pga3d, t: f32) -> pga3d {
  return add(a, mul(sub(b, a), t))
}

// slerp :: proc(a: pga3d, b: pga3d, t: f32) -> pga3d {
//   omega := math.acos( dot(a,b) )
//   return add(mul(a, math.sin((1.0-t)*omega)/math.sin(omega)), mul(b, math.sin(t*omega)/math.sin(omega)) )
// }

nlerp :: proc(a: pga3d, b: pga3d, t: f32) -> pga3d {
  return normalized(lerp(a, b, t))
}

// Primitives
rotor :: proc(angle: f32, line: pga3d) -> pga3d {
  return add(math.cos(angle/2.0), mul(math.sin(angle/2.0), normalized(line)))
}

translator :: proc(dist: f32, line: pga3d) -> pga3d {
  return 1.0+ dist/mul(2.0, line)
}

// not sure if needed if only used for static e0→e021
@(private="file")
make :: proc($f: f32, $idx: u32) -> (res: pga3d) {
  res[idx] = f
  return
}

// e0   : pga3d = make(1.0, 1)
// e1   : pga3d = make(1.0, 2)
// e2   : pga3d = make(1.0, 3)
// e3   : pga3d = make(1.0, 4)

// e01  : pga3d = outer(e0,e1)
// e02  : pga3d = outer(e0,e2)
// e03  : pga3d = outer(e0,e3)
// e12  : pga3d = outer(e1,e2)
// e31  : pga3d = outer(e3,e1)
// e23  : pga3d = outer(e2,e3)

// e123 : pga3d = outer(e1, e2, e3)
// e032 : pga3d = outer(e0, e3, e2)
// e013 : pga3d = outer(e0, e1, e3)
// e021 : pga3d = outer(e0, e2, e1)

// const basis
e0    : pga3d : {0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}
e1    : pga3d : {0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}
e2    : pga3d : {0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}
e3    : pga3d : {0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}

e0e1   : pga3d : {0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}
e0e2   : pga3d : {0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0}
e0e3   : pga3d : {0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0}
e1e2   : pga3d : {0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0}
e3e1   : pga3d : {0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0}
e2e3   : pga3d : {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0}

e1e2e3  : pga3d : {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0}
e0e3e2  : pga3d : {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0}
e0e1e3  : pga3d : {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0}
e0e2e1  : pga3d : {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0}

/* Helpers */
plane :: proc( a, b, c, d: f32) -> pga3d {
  return add( mul(a,e1), mul(b,e2), mul(c,e3,), mul(d,e0) )
}

point :: proc( x, y, z: f32) -> pga3d {
  return add( e1e2e3, mul(x,e0e3e2), mul(y,e0e1e3), mul(z,e0e2e1) )
}

circle :: proc( t, radius : f32, line: pga3d) -> pga3d {
  return mul( rotor((t * 2.0 * math.PI), line), translator(radius, mul(e1, e0)) )
}

torus :: proc( s, t, r1 : f32, l1 : pga3d, r2: f32, l2 : pga3d) -> pga3d {
  return mul( circle( s, r2, l2), circle( t, r1, l1) )
}

point_on_torus :: proc( s, t : f32) -> pga3d {
  return torus( s, t, 0.25, mul(e1, e2), 0.6, mul(e1, e3) )
}

/* Tests */
@test
test_pga3d :: proc(t: ^testing.T) {
  // todo! actually do some testing

  rot := rotor( math.PI / 2.0, outer( e1, e2 ) )

  ax_z := outer( e1, e2 )
  orig := outer( ax_z, e3)

  px := point( 1.0, 0.0, 0.0 )
  line := regress_pga3d( orig, px)

  p := plane(2,0,1,-3)

  rotated_plane := mul( rot, p, reverse(rot) )
  rotated_line := mul( rot, line, reverse(rot) )
  rotated_point := mul( rot, px, reverse(rot) )

  point_on_plane := mul( inner(p, px), p )

  fmt.println( "point:         ", px )
  fmt.println( "line:          ", rot )
  fmt.println( "plane:         ", p )
  fmt.println( "rotor:         ", rot )
  fmt.println( "rotated line:  ", rotated_line )
  fmt.println( "rotated point: ", rotated_point )
  fmt.println( "rotated plane: ", rotated_plane )
  fmt.println( "point on plane:", normalized(point_on_plane) )
  fmt.println( "point on torus:", point_on_torus(1.0, 0.0) )
}

@test
test_pga3d_fmt :: proc(t: ^testing.T) {
  p := plane(2,0,1,-3)

  fmt.printf( "%s: %v\n", "%v", p )
  fmt.printf( "%s: %f\n", "%f", p )
  fmt.printf( "%s: %T\n", "%T", p )
}