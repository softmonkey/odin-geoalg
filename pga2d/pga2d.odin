package bivector_pga2d

import "base:intrinsics"

import "core:math"
import "core:strings"
import "core:fmt"
import "core:reflect"
import "core:io"
import "core:testing"

pga2d     :: distinct [8]f32
pga2d_f16 :: distinct [8]f16
pga2d_f32 :: distinct [8]f32
pga2d_f64 :: distinct [8]f64

// Proc groups
add       :: proc { add_pga2d, add_multi_pga2d, adds_pga2d, sadd_pga2d, }
sub       :: proc { sub_pga2d, sub_multi_pga2d, ssub_pga2d, subs_pga2d, }
outer     :: proc { outer_pga2d, outer_multi_pga2d, }
mul       :: proc { mul_pga2d, mul_multi_pga2d, smul_pga2d, muls_pga2d, }
reverse   :: proc { reverse_pga2d, }
dual      :: proc { dual_pga2d, }
conj      :: proc { conjugate_pga2d, }
involute  :: proc { involute_pga2d, }
regress   :: proc { regress_pga2d, }
inner     :: proc { inner_pga2d, }


// fmt support
COLORED_FMT := true

@(private="file")
pga2d_basis_labels : [8]string = { "1","e0","e1","e2","e01","e20","e12","e012" }

@(init, private="file")
init :: proc() {
  fmt.set_user_formatters(&pga2d_formatters)
  err := fmt.register_user_formatter(type_info_of(pga2d).id, User_Formatter)
  assert(err == .None)
}

@(private="file")
pga2d_formatters: map[typeid]fmt.User_Formatter

@(private="file")
User_Formatter :: proc(fi: ^fmt.Info, arg: any, verb: rune) -> bool {
  m := cast(^pga2d)arg.data
  switch verb {
    case 'v':
      io.write_string(fi.writer, fmt_pga2d(m))
    case 'f':
      fmt.fmt_array(fi, m, 8, size_of(f32), type_info_of(f32),'f')
    case:
      return false
    }
    return true
}

@(private="file")
fmt_pga2d :: proc(a: ^pga2d) -> string {
  builder, err := strings.builder_make_len_cap(0,256,context.temp_allocator)
  if err != nil { return "" }
  defer strings.builder_destroy(&builder)

  tmpl := COLORED_FMT ? "%s\x1B[38;5;81m%.3f\x1B[38;5;214m%s\033[0m" : "%s%.3f%s"

  for val,idx in a {
    if math.abs(val) >= math.F32_EPSILON {
      fmt.sbprintf( &builder, tmpl, 
        len(builder.buf)>0 ? ", " : "", 
        val, 
        pga2d_basis_labels[idx] )
    }
  }
  if len(builder.buf) == 0 { return "\x1B[38;5;81m0.000\033[0m" }
  return strings.to_string(builder)
}

// Operations

// ~a
reverse_pga2d :: proc(a: pga2d) -> (res: pga2d) {
  res[0] =  a[0]
  res[1] =  a[1]
  res[2] =  a[2]
  res[3] =  a[3]
  res[4] = -a[4]
  res[5] = -a[5]
  res[6] = -a[6]
  res[7] = -a[7]
  return
}

// !a
dual_pga2d :: proc(a: pga2d) -> pga2d {
  return swizzle(a, 7,6,5,4,3,2,1,0)
}

conjugate_pga2d :: proc(a: pga2d) -> (res: pga2d) {
  res[0] =  a[0]
  res[1] = -a[1]
  res[2] = -a[2]
  res[3] = -a[3]
  res[4] = -a[4]
  res[5] = -a[5]
  res[6] = -a[6]
  res[7] =  a[7]
  return
}

involute_pga2d :: proc(a: pga2d) -> (res: pga2d){
  res[0] =  a[0]
  res[1] = -a[1]
  res[2] = -a[2]
  res[3] = -a[3]
  res[4] =  a[4]
  res[5] =  a[5]
  res[6] =  a[6]
  res[7] = -a[7]
  return
}

// *
mul_pga2d :: proc(a: pga2d, b: pga2d) -> (res: pga2d) {
  res[0] = b[0]*a[0]+b[2]*a[2]+b[3]*a[3]-b[6]*a[6]
  res[1] = b[1]*a[0]+b[0]*a[1]-b[4]*a[2]+b[5]*a[3]+b[2]*a[4]-b[3]*a[5]-b[7]*a[6]-b[6]*a[7]
  res[2] = b[2]*a[0]+b[0]*a[2]-b[6]*a[3]+b[3]*a[6]
  res[3] = b[3]*a[0]+b[6]*a[2]+b[0]*a[3]-b[2]*a[6]
  res[4] = b[4]*a[0]+b[2]*a[1]-b[1]*a[2]+b[7]*a[3]+b[0]*a[4]+b[6]*a[5]-b[5]*a[6]+b[3]*a[7]
  res[5] = b[5]*a[0]-b[3]*a[1]+b[7]*a[2]+b[1]*a[3]-b[6]*a[4]+b[0]*a[5]+b[4]*a[6]+b[2]*a[7]
  res[6] = b[6]*a[0]+b[3]*a[2]-b[2]*a[3]+b[0]*a[6]
  res[7] = b[7]*a[0]+b[6]*a[1]+b[5]*a[2]+b[4]*a[3]+b[3]*a[4]+b[2]*a[5]+b[1]*a[6]+b[0]*a[7]
  return
}
mul_multi_pga2d :: proc(a: ..pga2d) -> (res: pga2d) { 
  res = a[0]
  for i in a[1:] { 
    res = mul_pga2d(res, i)
  }
  return
}
smul_pga2d :: proc(a: f32, b: pga2d) -> pga2d { return a*b }
muls_pga2d :: proc(a: pga2d, b: f32) -> pga2d { return a*b }

// ^
outer_pga2d :: #force_inline proc(a: pga2d, b: pga2d) -> (res: pga2d) {
  res[0] = b[0]*a[0]
  res[1] = b[1]*a[0]+b[0]*a[1]
  res[2] = b[2]*a[0]+b[0]*a[2]
  res[3] = b[3]*a[0]+b[0]*a[3]
  res[4] = b[4]*a[0]+b[2]*a[1]-b[1]*a[2]+b[0]*a[4]
  res[5] = b[5]*a[0]-b[3]*a[1]+b[1]*a[3]+b[0]*a[5]
  res[6] = b[6]*a[0]+b[3]*a[2]-b[2]*a[3]+b[0]*a[6]
  res[7] = b[7]*a[0]+b[6]*a[1]+b[5]*a[2]+b[4]*a[3]+b[3]*a[4]+b[2]*a[5]+b[1]*a[6]+b[0]*a[7]
  return
}

outer_multi_pga2d :: proc(a: ..pga2d) -> (res: pga2d) { 
  res = a[0]
  for i in a[1:] { 
    res = outer_pga2d(res, i)
  }
  return
}

// &
regress_pga2d :: #force_inline proc(a: pga2d, b: pga2d) -> (res: pga2d) {
  res[7] = a[7]*b[7]
  res[6] = a[6]*b[7]+a[7]*b[6]
  res[5] = a[5]*b[7]+a[7]*b[5]
  res[4] = a[4]*b[7]+a[7]*b[4]
  res[3] = a[3]*b[7]+a[5]*b[6]-a[6]*b[5]+a[7]*b[3]
  res[2] = a[2]*b[7]-a[4]*b[6]+a[6]*b[4]+a[7]*b[2]
  res[1] = a[1]*b[7]+a[4]*b[5]-a[5]*b[4]+a[7]*b[1]
  res[0] = a[0]*b[7]+a[1]*b[6]+a[2]*b[5]+a[3]*b[4]+a[4]*b[3]+a[5]*b[2]+a[6]*b[1]+a[7]*b[0]
  return
}

// |
inner_pga2d :: #force_inline proc(a: pga2d, b: pga2d) -> (res: pga2d) {
  res[0] = b[0]*a[0]+b[2]*a[2]+b[3]*a[3]-b[6]*a[6]
  res[1] = b[1]*a[0]+b[0]*a[1]-b[4]*a[2]+b[5]*a[3]+b[2]*a[4]-b[3]*a[5]-b[7]*a[6]-b[6]*a[7]
  res[2] = b[2]*a[0]+b[0]*a[2]-b[6]*a[3]+b[3]*a[6]
  res[3] = b[3]*a[0]+b[6]*a[2]+b[0]*a[3]-b[2]*a[6]
  res[4] = b[4]*a[0]+b[7]*a[3]+b[0]*a[4]+b[3]*a[7]
  res[5] = b[5]*a[0]+b[7]*a[2]+b[0]*a[5]+b[2]*a[7]
  res[6] = b[6]*a[0]+b[0]*a[6]
  res[7] = b[7]*a[0]+b[0]*a[7]
  return 
}

// +
add_pga2d :: proc(a: pga2d, b: pga2d) -> pga2d { return a + b }
add_multi_pga2d :: proc(a: ..pga2d) -> (res: pga2d) { 
  res = a[0]
  for i in a[1:] { 
    res = add_pga2d(res, i)
  }
  return
}
sadd_pga2d :: proc(a: f32, b: pga2d) -> (res: pga2d) { 
  res = b
  res[0] += a
  return
}
adds_pga2d :: proc(a: pga2d, b: f32) -> (res: pga2d) { 
  res = a
  res[0] += b
  return
}

// -
sub_pga2d :: proc(a: pga2d, b: pga2d) -> pga2d { return a - b }
sub_multi_pga2d :: proc(a: ..pga2d) -> (res: pga2d) { 
  res = a[0]
  for i in a[1:] { 
    res = sub_pga2d(res, i)
  }
  return
}
ssub_pga2d :: proc(a: f32, b: pga2d) -> (res: pga2d) { 
  res = -b
  res[0] = a- b[0]
  return
}
subs_pga2d :: proc(a: pga2d, b: f32) -> (res: pga2d) { 
  res = -a
  res[0] = a[0] - b
  return
}

norm :: #force_inline proc(a: pga2d) -> f32 {
  return math.sqrt(math.abs(mul(a,conj(a))[0]))
}

inorm :: #force_inline proc(a: pga2d) -> f32 {
  return norm( dual(a) )
}

normalized :: #force_inline proc(a: pga2d) -> pga2d {
  return mul(a, 1.0/norm(a))
}

e0    : pga2d : {0,1,0,0,0,0,0,0}
e1    : pga2d : {0,0,1,0,0,0,0,0}
e2    : pga2d : {0,0,0,1,0,0,0,0}

e01   : pga2d : {0,0,0,0,1,0,0,0}
e20   : pga2d : {0,0,0,0,0,1,0,0}
e12   : pga2d : {0,0,0,0,0,0,1,0}

e012  : pga2d : {0,0,0,0,0,0,0,1}

@test
test_pga2d :: proc(t: ^testing.T) {
  fmt.println( "e0*e0:  ", mul(e0, e0) )
  fmt.println( "pss:    ", e012 )
  fmt.println( "pss*pss:", mul(e012, e012) )
}