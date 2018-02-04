function canvasApp() {
  $('#num').css({ color: 'black' })

  var dt = $("#in-res").val()
  var lines = dt.split('\n')
  var index = 0, lastLen = ''
  function doDraw() {
    if (index >= lines.length) {
      $('#num').css({ color: 'red' })
      return
    }

    var strs = lines[index].split(' ')

    if (strs[0] === lastLen) {
      // 如果和上次的结果相同，就跳过这一代
      ++index;
      setTimeout(doDraw, 0)
      return
    }

    ctx.clearRect(0, 0, 800, 800)

    $('#num').html(strs[0])
    $('#gen').html(index)

    var x, y;
    [x, y] = getPointFromDict(strs[1]);
    ctx.beginPath()
    ctx.moveTo(x, y);
    ctx.fillStyle = "rgb(0,255,0)"
    ctx.fillRect(x, y, 5, 5)
    ctx.fillStyle = "rgb(255,0,0)"
    for (var i = 2; i < strs.length; ++i) {
      [x, y] = getPointFromDict(strs[i]);
      ctx.fillRect(x, y, 5, 5)
      ctx.lineTo(x, y);
    }
    ctx.fillStyle = "rgb(0,0,255)"
    ctx.fillRect(x, y, 5, 5)
    ctx.stroke()

    ++index;
    lastLen = strs[0]

    setTimeout(doDraw, 500)
  }
  doDraw()
}

function translatePoint(x, y) {
  var res = [x * rate + zeroPosition[0], y * -1 * rate + zeroPosition[1]]
  //console.log(res)
  return res
}

function getPointFromDict(idx) {
  var item = points[idx]
  return translatePoint(item[0], item[1])
}

var max_DiffZero2 = 0
var points = {}
var zeroPosition = [400, 400]
var rate = 1

function initPoints() {
  for (var i in points) {
    var x = points[i][0]
    var y = points[i][1]
    var diffZero2 = x * x + y * y;
    max_DiffZero2 = Math.max(max_DiffZero2, diffZero2)
  }

  /*
  逻辑宽度/区域宽度 = 点的位置 / 点的显示位置
  点的显示位置 = 点的位置 * 区域宽度 / 逻辑宽度
  */
  rate = 400 / Math.sqrt(max_DiffZero2)
}

function convertInput() {
  var raw_input = $('#in-data').val();
  var in_data = raw_input.match(/\(.*?\)/gm)
  var res = ''
  var idx = 1
  for (var item of in_data) {
    // "(30,30)" -> "30,30"
    var s = item.substr(1, item.length - 2)
    var x, y;
    [x, y] = s.split(',')
    x = Number.parseFloat(x)
    y = Number.parseFloat(y)
    points[idx] = [x, y]
    res += idx + ' ' + x + ' ' + y + '\n'
    ++idx;
  }
  $('#out-data').val(in_data.length + '\n' + res)

  initPoints()
  ctx.fillStyle = "rgb(255,0,0)"
  for (i in points) {
    var x, y;
    [x, y] = translatePoint(points[i][0], points[i][1]);
    ctx.fillRect(x, y, 5, 5)
  }
}
