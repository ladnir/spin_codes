# otgen_patch(src dst anchor replacement): copy src to dst with the hook header
# prepended and `anchor` replaced. Fails unless the anchor occurs exactly once,
# so a change in the SPIN sources cannot silently drop a hook.
function(otgen_patch src dst anchor replacement)
  set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS ${src})
  file(READ ${src} text)
  string(FIND "${text}" "${anchor}" first)
  string(FIND "${text}" "${anchor}" last REVERSE)
  if(first EQUAL -1 OR NOT first EQUAL last)
    message(FATAL_ERROR "otgen_patch: anchor must occur exactly once in ${src}: ${anchor}")
  endif()
  string(REPLACE "${anchor}" "${replacement}" text "${text}")
  file(WRITE ${dst}.tmp "#include \"otgen_hook.h\"\n${text}")
  configure_file(${dst}.tmp ${dst} COPYONLY)
endfunction()
