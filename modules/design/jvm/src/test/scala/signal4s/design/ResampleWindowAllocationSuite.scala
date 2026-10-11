package signal4s.design

import gale.linalg.Vec

class ResampleWindowAllocationSuite extends munit.FunSuite:
  private def ok[E,A](v: Either[E,A]): A=v.fold(e=>fail(e.toString),identity)
  test("ordinary512-output windows avoid per-output arbitrary-precision clock allocation"):
    val plan=ok(ResamplePolyPlan(2,3));val window=ok(plan.inputWindow(6007L,1000L,512))
    val input=Vec.tabulate(window.count)(i=>(i%31-15)*1e-5)
    for _<-0 until 5 do {val _=ok(plan.processWindow(input,window.start,6007L,1000L,512));()}
    val bean=java.lang.management.ManagementFactory.getThreadMXBean.asInstanceOf[com.sun.management.ThreadMXBean]
    assert(bean.isThreadAllocatedMemorySupported);bean.setThreadAllocatedMemoryEnabled(true)
    val thread=Thread.currentThread().threadId();val before=bean.getThreadAllocatedBytes(thread)
    val result=ok(plan.processWindow(input,window.start,6007L,1000L,512))
    val allocated=bean.getThreadAllocatedBytes(thread)-before
    assert(allocated<131072L,s"one ordinary512-output window allocated $allocated bytes")
    assertEquals(result.length,512)
