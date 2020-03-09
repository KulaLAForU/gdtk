#! /usr/bin/env ruby
# mms-euler-test.rb
# Run the method of manufactured solutions test case to exercise
# the Euler equations using the unstructured (steady-state) solver.
# Our tests involve checking the density norms against
# known good values. These known good values come from an order
# of accuracy test.
# KAD, 2020-02-29;
#
require 'test/unit'
require 'open3'

class TestMMS_Euler < Test::Unit::TestCase
  def test_0_prep
    cmd = "cp case-euler.txt case.txt"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    cmd = "python3 make_lua_files.py"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    cmd = "e4shared --prep --job=mms"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
  end

  def test_1_run
    cmd = "e4zsss --job=mms"
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    steps = 0
    lines = o.split("\n")
    lines.each do |txt|
      if txt.match('STEP=') then
        items = txt.split(' ')
        steps = items[1].to_i
      end
    end
    assert((steps - 113).abs < 3, "Failed to take correct number of steps.")
  end

  def test_2_check_norms
    cmd = 'e4shared --job=mms --post --tindx-plot=last --ref-soln=ref-soln.lua --norms="rho"'
    o, e, s = Open3.capture3(*cmd.split)
    assert_equal(s.success?, true)
    lines = o.split("\n")
    names = ['L1', 'L2', 'Linf']
    refs = {'L1'=>1.730145442812708149e-03, 'L2'=>2.436181746051479212e-03, 'Linf'=>6.587021426123662593e-03}
    values = {'L1'=>0.0, 'L2'=>0.0, 'Linf'=>0.0}
    lines.each do |txt|
      if txt.match('Linf')
        # Must be our data line.
        items = txt.strip().split(' ')
        values['L1'] = items[1].to_f
        values['L2'] = items[3].to_f
        values['Linf'] = items[5].to_f
        break
      end
    end
    names.each do |n|
      assert((values[n]-refs[n]).abs/refs[n] < 1.0e-3, "Failed to see correct #{n}.")
    end
  end
end